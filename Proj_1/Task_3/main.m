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
% SI-mm units used, rho does not matter, everything already inside

pbFC_name = 'my_pbFC';     % Points, bars, forces and boundary conditons

%% Initialize EDACFEM ****************************
pbFC = load(pbFC_name).pbFC;  %Load points, bars, forces and boundary conditons

% Parameters for solver
switch_importMethod = 'pbFC'; %'script_full', 'script_simplified', 'pbFC'
switch_outputMode = 'silent'; %'verbose', 'silent'
% Creating FE model, correct values in pbFC already
[fem,opts] = import_model(pwd,switch_importMethod,switch_outputMode,pbFC);

% Compute structure response and visualize
[fem] = performFEM(fem,opts);



%% Compute sensitivities *************************


x_matrix=[-500 -500 ; 500 500; -500 500; 500 -500; 0 0]'; % x1 to x5

f_matrix = zeros(1,5);
for i=1:5
    f_matrix(i) = Objective(x_matrix(:,i)',fem,opts);
end

[max_error,i_max_error]=max(abs(f_matrix'-single(f_matrix')));
dx = double(sqrt(max_error / abs(Objective(x_matrix(:,i_max_error)',fem,opts))));



x0 = [200, 100];

s = sensitivity(x0,dx,fem,opts);
format long e
disp('Sensitivitys for x0 = [200, 100]:')
disp(s)


%% Compute Sensitivity with first order FD
function sens = sensitivity(x,dx,fem,opts)
    % value at point
    f_0 = Objective(x,fem,opts);

    fx = (Objective(x+[dx,0],fem,opts)-f_0)/dx;
    fy = (Objective(x+[0,dx],fem,opts)-f_0)/dx;

    sens = [fx,fy];

end


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


