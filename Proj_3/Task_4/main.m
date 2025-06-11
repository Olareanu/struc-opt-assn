%***********************************
%
%   Structural Optimization 2025
%   Olareanu Alexandru
%   Proj 3 task 4
%
%% ***********************************************
clear all;  %clear workspace
close all;  %close figures
clc;
%% ***********************************************

% Input Parameters
h0 = 0.5;
yf = 1;

% Measure execution time
tic;  % Start timer
[V_optimal, A_optimal] = optimizeTruss(h0, yf);
elapsed_time = toc;  % Stop timer

% Print the elapsed time
fprintf('Time taken to run optimizeTruss: %.2f seconds\n', elapsed_time);





%% Truss Generation and Optimization

function [fval, A_optimal] = optimizeTruss(h0, yf)

% Mesh generation (as per Exercise) *************
drectangle = @(p,x1,x2,y1,y2) -min(min(min(-y1+p(:,2),y2-p(:,2)),-x1+p(:,1)),x2-p(:,1));
fd = @(p)  drectangle(p,0,2,0,2);           %Domain
fh = @(p) ones(size(p,1),1);                % Uniform mesh
bounding_box=[0,0;2,2];                     %Domain bounding-box
p_keep=[0,0;0,2;2,0;2,2;0,0.5;0,1.5;2,yf];  %Points that are always kept in the generated mesh
[p,t,b,L] = distmeshSO( fd, fh, h0, bounding_box, p_keep );
patch( 'vertices', p, 'faces', t, 'facecolor', [.9, .9, .9] )
title('Initial Mesh');
axis equal


% Boundary conditions ****************************
eps_=1e-5; %% Geometric search precision as distmesh is numeric
% Supports **************************************
locsup = find(abs(p(:,1))<=eps_ & (abs(p(:,2)-0.5)<=eps_ | abs(p(:,2)-1.5)<=eps_));
assert(length(locsup) == 2, 'Expected exactly 2 support points but found %d', length(locsup));
locsup=[locsup*2;locsup*2-1]; % x and y DOF locked
% Loads *****************************************
locf = find(abs(p(:,1)-2)<=eps_ & abs(p(:,2)-yf)<=eps_); %find point for load at (2,yf)
assert(length(locf) == 1, 'Expected exactly 1 load point but found %d', length(locf));
locf=[locf*2;locf*2-1]; % x DOF loaded


% Optimization using fmincon *********************
% Initial design
A_initial = 10*ones(size(b,1),1);  % Initial cross-section areas

% Create objective function (minimizing volume)
objectiveFunction = @(A) calculateVolume(L, A);

% Options for fmincon
options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'sqp', 'MaxFunctionEvaluations', 1e5, ...
    'SpecifyObjectiveGradient', true);

% Create nonlinear constraint function for stress limits
nonlcon = @(A) stressConstraints(p, t, b, L, locf, locsup, A);

% No linear constraints
A_ineq = [];
b_ineq = [];
A_eq = [];
b_eq = [];

% For now, just use wide bounds to ensure positive areas
lb = 1e-6*ones(size(A_initial));  % Lower bound
ub = 20*ones(size(A_initial));   % Upper bound

% Run optimization
[A_optimal, fval] = fmincon(objectiveFunction, A_initial, A_ineq, b_ineq, A_eq, b_eq, lb, ub, nonlcon, options);

fprintf('\nOptimization Results:\n');
fprintf('Minimum volume: %f\n', fval);

% Compute the optimized structure for both load cases
[p_deformed_opt1, Sigma_opt1, N_opt1, p_deformed_opt2, Sigma_opt2, N_opt2] = computeFe(p, t, b, L, locf, locsup, A_optimal);

% Print optimized cross-sectional areas and stresses for both load cases
fprintf('\nOptimized cross-sectional areas and stresses for both load cases:\n');
member_indices = (1:length(Sigma_opt1))';
table(member_indices, A_optimal, N_opt1, Sigma_opt1, N_opt2, Sigma_opt2, ...
    'VariableNames', {'Member', 'Area', 'Force_LC1', 'Stress_LC1', 'Force_LC2', 'Stress_LC2'})

% Find the node indices for visualization (convert from DOFs back to node numbers)
locsup_nodes = unique(ceil(locsup/2));
locf_node = unique(ceil(locf(1)/2));

% Plot initial structure
Plot_Structure(2, p, b, A_optimal, 'Optimized Structure (Undeformed) - Volume:', fval, 5, 0.1, locsup_nodes(1), locsup_nodes(2), locf_node);

% Plot optimized structure for Load Case 1
Plot_Structure(3, p_deformed_opt1, b, A_optimal, 'Optimized Structure (Deformed - Load Case 1) - Volume:', fval, 5, 0.1, locsup_nodes(1), locsup_nodes(2), locf_node);

% Plot optimized structure for Load Case 2
Plot_Structure(4, p_deformed_opt2, b, A_optimal, 'Optimized Structure (Deformed - Load Case 2) - Volume:', fval, 5, 0.1, locsup_nodes(1), locsup_nodes(2), locf_node);

end


%% Stress Constraint Function *********************
function [c, ceq] = stressConstraints(p, t, b, L, locf, locsup, A)
% Calculate stresses for the current design
[~, Sigma_1, ~, ~, Sigma_2, ~] = computeFe(p, t, b, L, locf, locsup, A);

Sigma_y = 1;  % Yield stress

% Inequality constraints: stress must be less than yield stress for both load cases
c1 = abs(Sigma_1) - Sigma_y;
c2 = abs(Sigma_2) - Sigma_y;
c = [c1; c2]; % Combine constraints from both load cases

% No equality constraints
ceq = [];
end


%% Unified FE Function for Both Load Cases *************************************
function [p_deformed_1, Sigma_1, N_1, p_deformed_2, Sigma_2, N_2] = computeFe(p, t, b, L, locf, locsup, A)

% FE Parameters
E = 1;  % Young's modulus

% Load Case 1
Fx1 = 0;   % Load in x-direction
Fy1 = -2;  % Load in y-direction

% Load Case 2
Fx2 = 1;   % Load in x-direction
Fy2 = 0;   % Load in y-direction

% FE Assembly (done only once)
n = size(p,1)*2;              % Full DOFs
S = sparse(diag(E.*A./L));    % Local stiffness
B = sparse(B_generator(p, b));
Kg = B*S*B';                  % Full K
K = Kg;
K(locsup,:) = [];             % Remove constrained DOFs
K(:,locsup) = [];             % Remove constrained DOFs

% Setup load vectors for both cases
f1 = sparse(zeros(n,1));      % Load vector for case 1
f1(locf(1),:) = Fx1;          % x-component load
f1(locf(2),:) = Fy1;          % y-component load
f1(locsup,:) = [];            % Remove constrained DOFs

f2 = sparse(zeros(n,1));      % Load vector for case 2
f2(locf(1),:) = Fx2;          % x-component load
f2(locf(2),:) = Fy2;          % y-component load
f2(locsup,:) = [];            % Remove constrained DOFs

% Solve for both load cases
u1 = K\f1;
u2 = K\f2;

% Remove constrained DOFs from B for both cases
B_reduced = B;
B_reduced(locsup,:) = [];

% Compute member forces and stresses for both cases
dl1 = -(B_reduced'*u1);       % Member elongation case 1
N_1 = E.*A./L.*dl1;           % Member force case 1
Sigma_1 = N_1./A;             % Member stress case 1

dl2 = -(B_reduced'*u2);       % Member elongation case 2
N_2 = E.*A./L.*dl2;           % Member force case 2
Sigma_2 = N_2./A;             % Member stress case 2

% Merge displacements with supports for both cases
pr = 1:n;
pr(locsup) = [];

U1 = sparse(zeros(n,1));
U1(pr,1) = u1;
p_deformed_1 = [p(:,1)+U1(1:2:n), p(:,2)+U1(2:2:n)];

U2 = sparse(zeros(n,1));
U2(pr,1) = u2;
p_deformed_2 = [p(:,1)+U2(1:2:n), p(:,2)+U2(2:2:n)];

end


%% Volume as Objective Function *********************
function [V, gradV] = calculateVolume(L, A)
% Objective function
V = sum(L .* A);

% Gradient of the objective function
gradV = L;
end