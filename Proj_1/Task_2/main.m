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


%% Set up optimization ***************************
x0 = [0,0]; % Start point

% History arrays
global xHistory FHistory

xHistory = [];
FHistory = [];

f = @(x)Objective(x,fem,opts);

optimizer_options = optimoptions(@fminunc,'Display','iter-detailed','Algorithm','quasi-newton', 'PlotFcn',{@optimplotfval, @optimplotx});
[x_optimal,f_optimal] = fminunc(f,x0,optimizer_options);

disp('f_optimal:');
disp(f_optimal);
disp('x_optimal:');
disp(x_optimal);

% final X, Y position of node 4
final_positions = x_optimal; 
final_positions = [400, 100] + x_optimal;
disp('final_positions_node_4:');
disp(final_positions);


% compute resulting structure one last time and show it

% Set position of node 4 to that coordinates
fem.p(4,1:2)=final_positions;

% Update struct and compute
[fem,opts] = setParamsBeamFEM(fem,opts);
[fem] = performFEM(fem,opts);

% visualize
viz2D3D_line_deformed(fem,opts,1,'Mag'); 
viz2D3D_line_stresses(fem,opts,1);

% Get the reaction forces
F = fem.sol{1,1}.F;

R1 = F(1,1:3);
R2 = F(2,1:3);
R3 = F(3,1:3);

% Compute objective function
f = sqrt(norm(R1)^2 + norm(R2)^2 + norm(R3)^2);
% disp('f_final:');
% disp(f);


% Plot History
figure;
subplot(2,1,1);
plot(1:length(FHistory), FHistory, 'o-');
title('Objective Function');
xlabel('iter');
ylabel('f(x)');

subplot(2,1,2);
plot(1:length(xHistory), xHistory(:,1), 'o-', 'DisplayName', 'x_1');
hold on;
plot(1:length(xHistory), xHistory(:,2), 's-', 'DisplayName', 'x_2');
title('x Values');
xlabel('iter');
ylabel('x Values');
legend('show');


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

% Append current x and f to history
global xHistory FHistory
xHistory = [xHistory; x(:)'];
FHistory = [FHistory; f];

end


