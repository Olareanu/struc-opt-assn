%************************************************
%
%   Structural Optimization 2025
%   Olareanu Alexandru
%   Proj 2 task 4
%
%% **********************************************
clear all;  %clear workspace
close all;  %close figures
clc;
%% Input parameters *****************************
% Dimensions follow the naming convention of the diagram

% Load pre-computed data for surrogate models
load results.mat
load weightsDisp.mat
load weightsMises.mat

% Define the objective function (panel volume to be minimized)
fun = @(x)panel_volume(x(1),x(2));

% Set optimization parameters
x0 = [5, 50];        % Initial values for [panel thickness, stiffener height]
lowerBound = [5, 50];        % Lower bounds for [tp, hs]
upperBound = [20, 120];      % Upper bounds for [tp, hs]

% Define constraints for displacement and stress
nonlcon = @(x)constraints(x(1),x(2),weigthsDisplacement,weigthsMises);

% Configure optimization options
options = optimoptions(@fmincon,...
                      'Algorithm', 'interior-point',...
                      'Display', 'iter-detailed',...
                      'MaxIterations', 1500,...
                      'PlotFcn', 'optimplotfval');

% Run optimization using fmincon
[x, fval] = fmincon(fun, x0, [], [], [], [], lowerBound, upperBound, nonlcon, options);

% Extract optimal design parameters
tp_opt = x(1);       % Optimal panel thickness (mm)
hs_opt = x(2);       % Optimal stiffener height (mm)
V_opt = fval;        % Resulting minimum volume (mm^3)

% Evaluate performance at optimal design point
% Using surrogate models
maxDisplacementSurr = surrogate(tp_opt, hs_opt, weigthsDisplacement);
maxVonMisesSurr = surrogate(tp_opt, hs_opt, weigthsMises);

% Verify with FE Model
[maxDisplacementFE, maxVonMisesFE, ~] = feComputation(tp_opt, hs_opt, true);

% Display optimization results
disp("Optimum found at t_p = " + num2str(tp_opt) + " mm and h_s = " + num2str(hs_opt) + " mm with V = " + num2str(fval) + "mm^3")
disp("Maximum displacement: FE: " + num2str(maxDisplacementFE) + " mm, Surrogate: " + num2str(maxDisplacementSurr) + " mm")
disp("Maximum von Mises: FE: " + num2str(maxVonMisesFE) + " MPa, Surrogate: " + num2str(maxVonMisesSurr) + " MPa")


%% **********************************************

% compute value from surrogate model
function value = surrogate(tp,hs,weights)
        value = [tp,hs,tp^2,hs^2,tp*hs,1] * weights;
    end


function [const, snedOver] = constraints(tp,hs,weightsDisplacement,weightsMises)
    const = [surrogate(tp,hs,weightsDisplacement)-3;
        surrogate(tp,hs,weightsMises)-150];
    snedOver = [];
    end

% compute pannel volume
function pannelVolume = panel_volume(tp,hs)

    L = 1000;
    W = 1000;
    bs = 5;

    pannelVolume = L*W*tp + 4*L*hs*bs;
    
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