%***********************************
%
%   Structural Optimization 2025
%   Olareanu Alexandru
%   Proj 1 task 1
%
%% ***********************************************
clear all;  %clear workspace
close all;  %close figures
clc;
%% Input parameters ******************************
% SI-mm units used, rho does not matter

pbFC_name = 'my_pbFC';     % Points, bars, forces and boundary conditons
A0=15;                  % Cross-section
E0=210000;              % Young's modulus
ro0=1;                  % Density

%% Initialize EDACFEM ****************************
pbFC = load(pbFC_name).pbFC;  %Load points, bars, forces and boundary conditons

% Parameters for solver
switch_importMethod = 'pbFC'; %'script_full', 'script_simplified', 'pbFC'
switch_outputMode = 'silent'; %'verbose', 'silent'
% Creating FE model, correct values in pbFC already
[fem,opts] = import_model(pwd,switch_importMethod,switch_outputMode,pbFC);

% Compute structure response and visualize
[fem] = performFEM(fem,opts);
viz2D3D_line_deformed(fem,opts,1,'Mag')
disp('Structural response in initial config:')
disp('Forces matrix:')
disp(fem.sol{1,1}.F)
disp('Displacements matrix:')
disp(fem.sol{1,1}.u)


%% Objective function ****************************
x0 = [-200,-200];
[f] = Objective(x0, fem, opts);

disp('Objective function at x0 = -200, -200')
display(f);




%% Compute objective function
function f = Objective(x, fem, opts)

% Compute new positions of node 4 X and Y
x = [400,100] + x;

% Set position of node 4 to that coordinates
fem.p(4,1:2)=x;

% Optionally print new points matrix
% fem.p

% Update struct and compute
[fem,opts] = setParamsBeamFEM(fem,opts);
[fem] = performFEM(fem,opts);

% optionally visualize
% viz2D3D_line_deformed(fem,opts,1,'Mag')

% Get the reaction forces
F = fem.sol{1,1}.F;

R1 = F(1,1:3);
R2 = F(2,1:3);
R3 = F(3,1:3);

% Compute objective function
f = sqrt(norm(R1)^2 + norm(R2)^2 + norm(R3)^2);


end


