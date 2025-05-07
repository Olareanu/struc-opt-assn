%************************************************
%
%   Structural Optimization 2025
%   Olareanu Alexandru
%   Proj 2 task 3
%
%% **********************************************
clear all;  %clear workspace
close all;  %close figures
clc;
%% Input parameters *****************************
% Dimensions follow the naming convention of the diagram

load results.mat
n = size(results,1);

% init empty matrixes and vectors
featureMatrix = zeros(n,6);
feDisp = zeros(n,1);
feMises = zeros(n,1);

% fill matrixes and vectors with data from table
for i = 1:n
    tp = results.tp{i};
    hs = results.hs{i};
    featureMatrix(i,:) = [tp,hs,tp^2,hs^2,tp*hs,1];
    feDisp(i) = results.maxDisplacement{i};
    feMises(i) = results.maxVonMises{i};
end

weigthsDisp = featureMatrix\feDisp;
weigthsMises = featureMatrix\feMises;
% weigthsDisplacement = lsqr(featureMatrix,FE_displacement);
% weigthsMises = lsqr(featureMatrix,FE_mises,1e-8,100);


save weightsDisplacement.mat weigthsDisp
save weightsMises.mat weigthsMises


% Original Response:
tpOriginal = 10;
hsOriginal = 90;
[maxDisplacementFE,maxVonMisesFE,~] = feComputation(tpOriginal,hsOriginal,true);

% Surrogate Model:
maxDisplacementSurrogate = surrogate(tpOriginal,hsOriginal,weigthsDisp);
maxVonMisesSurrogate = surrogate(tpOriginal,hsOriginal,weigthsMises);


disp("Max disp: FE: " + num2str(maxDisplacementFE) + " mm, surr: " + num2str(maxDisplacementSurrogate) + " mm")
disp("Max mises: FE: " + num2str(maxVonMisesFE) + " MPa, surr: " + num2str(maxVonMisesSurrogate) + " MPa")


%% **********************************************

% compute value from surrogate model
function value = surrogate(tp,hs,weights)
        value = [tp,hs,tp^2,hs^2,tp*hs,1] * weights;
    end


%% **********************************************

function [maxDisplacement,maxVonMises,Results] = feComputation(tp,hs,plot)
    
    L = 1000;
    W = 1000;
    bs = 5;

    % Inspired by Exercise 9
    structuralmodel = femodel("AnalysisType","structuralStatic");

    % Points for sketch, tedious
    points = [0,0;
         W/5-bs/2,0;
         W/5-bs/2,hs;
         W/5+bs/2,hs;
         W/5+bs/2,0;
         2*W/5-bs/2,0;
         2*W/5-bs/2,hs;
         2*W/5+bs/2,hs;
         2*W/5+bs/2,0;
         3*W/5-bs/2,0;
         3*W/5-bs/2,hs;
         3*W/5+bs/2,hs;
         3*W/5+bs/2,0;
         4*W/5-bs/2,0;
         4*W/5-bs/2,hs;
         4*W/5+bs/2,hs;
         4*W/5+bs/2,0;
         W,0;
         W,-tp;
         0,-tp];

    noPoints = size(points,1);

    Beam = [3;
             noPoints;
             points(:,1);
             points(:,2)];

    % Creating sketch out of points
    gm = [Beam];
    sf = '(Beam)';
    ns = char('Beam');
    ns = ns';
    [g,t] = decsg(gm,sf,ns);
    [g,t] = csgdel(g,t); 
    gm=geometryFromEdges(g);


    % Visualize sketch
    if plot
        figure
        pdegplot(g,'VertexLabels','on','EdgeLabels','on');
        title('2D Geometry')
    end

    % Extrude the sketch into a 3D object
    g3D = extrude(gm,L);

    % Visualize the 3D Gemetry
    if plot
        figure
        pdegplot(g3D,'faceLabels','on','FaceAlpha',0.25,'EdgeLabels','on','VertexLabels','on','CellLabels','on');
        title('3D Geometry');
    end

    % FE sim: meshing
    structuralmodel.Geometry = g3D;  %Assign geometry
    structuralmodel=structuralmodel.generateMesh('GeometricOrder','quadratic','Hmax',20,'Hmin',2.5);
    if plot
        figure
        pdemesh(structuralmodel,'FaceAlpha',0.5); 
        title('FE Mesh');
    end

    % FE sim: material properties
    E = 210000; % N/mm^2 MPa
    nu  = 0.3;
    structuralmodel.MaterialProperties(1)=materialProperties('YoungsModulus',E,'PoissonsRatio',nu);

    % FE sim: BCs
    % Edges to be constrained: E60, E19, E39, E59
    structuralmodel.EdgeBC([19 39 59 60]) = edgeBC("Constraint","fixed");

    % FE sim: Loads
    % Loaded face: F21, with p = rho * g * h
    p = 1000 * 10 * 10; % Pa
    p = p * 1e-6; % MPa
    structuralmodel.FaceLoad(21) = faceLoad("Pressure",p);

    % FE sim: solving
    % structuralmodel.solve();
    Results = solve(structuralmodel); %Get the results
    if plot
        figure
        scale=20; %Scaling factor for displacements visualization
        pdeplot3D(Results.Mesh,'ColorMapData',Results.Displacement.y, ...
                                  'Deformation',Results.Displacement, ...
                                  'DeformationScaleFactor',scale);
        title('Displacement z-direction');
        figure
        pdeplot3D(Results.Mesh,'ColorMapData',Results.VonMisesStress,'Deformation',Results.Displacement,'DeformationScaleFactor',scale);
        title('Von Mises Stress');
    end

    
    % Return maximums
    maxDisplacement = max(Results.Displacement.y);
    maxVonMises = max(Results.VonMisesStress);
end