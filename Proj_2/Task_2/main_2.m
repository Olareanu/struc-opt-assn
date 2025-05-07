%************************************************
%
%   Structural Optimization 2025
%   Olareanu Alexandru
%   Proj 2 task 2
%
%% **********************************************
clear all;  %clear workspace
close all;  %close figures
clc;
%% Input parameters *****************************
% Dimensions follow the naming convention of the diagram

n = 5;

tpSweep = linspace(5,20,n);
hsSweep = linspace(50,120,n);

results = table;

for i = 1:n
    tp = tpSweep(i);
    for j = 1:n
        disp("Running sim nr " + num2str((n)*(i-1) + j) + " out of " + num2str(n^2))

        hs = hsSweep(j);
        [maxDisplacement, maxVonMises,~] = feComputation(tp,hs,false);
        results.tp{n*(i-1)+j} = tp;
        results.hs{n*(i-1)+j} = hs;
        
        results.maxDisplacement{n*(i-1)+j} = maxDisplacement;
        results.maxVonMises{n*(i-1)+j} = maxVonMises;
    end
end

save('results.mat', 'results')

% CSV format
writetable(results, 'results.csv')

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