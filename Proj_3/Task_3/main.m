%***********************************
%
%   Structural Optimization 2025
%   Olareanu Alexandru
%   Proj 3 task 3
%
%% ***********************************************
clear all;  %clear workspace
close all;  %close figures
clc;
%% ***********************************************

% Input Parameters
h0 = 0.5;
yf = 1;

%% Mesh generation (as per Exercise) *************
drectangle = @(p,x1,x2,y1,y2) -min(min(min(-y1+p(:,2),y2-p(:,2)),-x1+p(:,1)),x2-p(:,1));
fd = @(p)  drectangle(p,0,2,0,2);           %Domain
fh = @(p) ones(size(p,1),1);                % Uniform mesh
bounding_box=[0,0;2,2];                     %Domain bounding-box
p_keep=[0,0;0,2;2,0;2,2;0,0.5;0,1.5;2,yf];  %Points that are always kept in the generated mesh
[p,t,b,L] = distmeshSO( fd, fh, h0, bounding_box, p_keep );
patch( 'vertices', p, 'faces', t, 'facecolor', [.9, .9, .9] )
title('Initial Mesh');
axis equal


%% Boundary conditions ****************************
eps_=1e-5; %% Geometric search precision as distmesh is numeric
% Supports **************************************
locsup = find(abs(p(:,1))<=eps_ & (abs(p(:,2)-0.5)<=eps_ | abs(p(:,2)-1.5)<=eps_));
assert(length(locsup) == 2, 'Expected exactly 2 support points but found %d', length(locsup));
locsup=[locsup*2;locsup*2-1]; % x and y DOF locked
% Loads *****************************************
locf = find(abs(p(:,1)-2)<=eps_ & abs(p(:,2)-yf)<=eps_); %find point for load at (2,yf)
assert(length(locf) == 1, 'Expected exactly 1 load point but found %d', length(locf));
locf=[locf*2;locf*2-1]; % x DOF loaded


%% Optimization using fmincon *********************
% Initial design
A_initial = 10*ones(size(b,1),1);  % Initial cross-section areas

% Create objective function (minimizing volume)
objectiveFunction = @(A) calculateVolume(L, A);

% Options for fmincon
options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'sqp');

% Create nonlinear constraint function for stress limits
nonlcon = @(A) stressConstraints(p, t, b, L, locf, locsup, A);

% No linear constraints
A_ineq = [];
b_ineq = [];
A_eq = [];
b_eq = [];

% For now, just use wide bounds to ensure positive areas
lb = 1e-6*ones(size(A_initial));  % Lower bound
ub = 100*ones(size(A_initial));   % Upper bound

% Run optimization
[A_optimal, fval] = fmincon(objectiveFunction, A_initial, A_ineq, b_ineq, A_eq, b_eq, lb, ub, nonlcon, options);

fprintf('\nOptimization Results:\n');
fprintf('Minimum volume: %f\n', fval);

% Compute the optimized structure
[p_deformed_opt, Sigma_opt, N_opt] = computeFe(p, t, b, L, locf, locsup, A_optimal);

% Print optimized cross-sectional areas and stresses
fprintf('\nOptimized cross-sectional areas and stresses:\n');
member_indices = (1:length(Sigma_opt))';
table(member_indices, A_optimal, Sigma_opt, N_opt, 'VariableNames', {'Member', 'Area', 'Stress', 'Force'})





% Find the node indices for visualization (convert from DOFs back to node numbers)
locsup_nodes = unique(ceil(locsup/2));
locf_node = unique(ceil(locf(1)/2));

% Plot initial structure
Plot_Structure(3, p, b, A_optimal, 'Optimized Structure (Undeformed) - Volume:', fval, 5, 0.1, locsup_nodes(1), locsup_nodes(2), locf_node);

% Plot optimized structure
Plot_Structure(4, p_deformed_opt, b, A_optimal, 'Optimized Structure (Deformed) - Volume:', fval, 5, 0.1, locsup_nodes(1), locsup_nodes(2), locf_node);




%% Fully Stressed Design to compare **************
fprintf('\n\n========== Fully Stressed Design ==========\n');

% Initial design
A_fsd = 10*ones(size(b,1),1);  % Initial cross-section areas
max_iter = 20;                 % Maximum iterations
tol = 1e-4;                    % Convergence tolerance
sigma_y = 1;                   % Yield stress

% Iteration history
A_history = zeros(length(A_fsd), max_iter+1);
A_history(:,1) = A_fsd;

% FSD iteration
for iter = 1:max_iter
    % Calculate stresses
    [p_deformed_fsd, Sigma_fsd, N_fsd] = computeFe(p, t, b, L, locf, locsup, A_fsd);
    
    % Calculate new areas
    A_new = abs(N_fsd) / sigma_y;
    
    % Update areas (with relaxation to improve convergence)
    relax = 0.5;  % Relaxation factor
    A_fsd = relax*A_new + (1-relax)*A_fsd;
    
    % Apply lower bound to prevent zero areas
    A_fsd = max(A_fsd, 1e-6);
    
    % Store history
    A_history(:,iter+1) = A_fsd;
    
    % Check convergence
    change = norm(A_new - A_fsd) / norm(A_fsd);
    volume = sum(L .* A_fsd);
    
    fprintf('Iteration %d: Volume = %.6f, Max change = %.6f\n', iter, volume, change);
    
    if change < tol
        fprintf('FSD converged after %d iterations\n', iter);
        break;
    end
end

% Final calculation with optimized areas
[p_deformed_fsd, Sigma_fsd, N_fsd] = computeFe(p, t, b, L, locf, locsup, A_fsd);
volume_fsd = sum(L .* A_fsd);

fprintf('\nFully Stressed Design Results:\n');
fprintf('Final volume: %f\n', volume_fsd);

% Print optimized cross-sectional areas and stresses
fprintf('\nFSD cross-sectional areas and stresses:\n');
member_indices = (1:length(Sigma_fsd))';
table(member_indices, A_fsd, Sigma_fsd, N_fsd, 'VariableNames', {'Member', 'Area', 'Stress', 'Force'})

% Plot FSD structure
Plot_Structure(5, p, b, A_fsd, 'FSD Structure (Undeformed) - Volume:', volume_fsd, 5, 0.1, locsup_nodes(1), locsup_nodes(2), locf_node);
Plot_Structure(6, p_deformed_fsd, b, A_fsd, 'FSD Structure (Deformed) - Volume:', volume_fsd, 5, 0.1, locsup_nodes(1), locsup_nodes(2), locf_node);






% Compare fmincon and FSD results
fprintf('\nComparison of optimization results:\n');
fprintf('fmincon volume: %.6f\n', fval);
fprintf('FSD volume: %.6f\n', volume_fsd);
fprintf('Difference: %.6f (%.2f%%)\n', volume_fsd - fval, 100*(volume_fsd - fval)/fval);













%% Stress Constraint Function *********************
function [c, ceq] = stressConstraints(p, t, b, L, locf, locsup, A)
% Calculate stresses for the current design
[~, Sigma, ~] = computeFe(p, t, b, L, locf, locsup, A);

Sigma_y = 1;  % Yield stress

c = abs(Sigma) - Sigma_y; % Inequality constraints: stress must be less than yield stress

% No equality constraints
ceq = [];
end


%% FE ********************************************
function [p_deformed, Sigma, N] = computeFe(p, t, b, L, locf, locsup, A)

% FE Parameters
F=.5;                       %Load
E=1;                        %Young's modulus


% FE Assembly and Solving
n=size(p,1)*2;              %Full DOFs
S=sparse(diag(E.*A./L));    %Local stiffness
B = sparse(B_generator(p, b));
Kg=B*S*B';                  %Full K
K=Kg;
K(locsup,:)=[];             %Remove constrained DOFs
K(:,locsup)=[];             %Remove constrained DOFs
f=sparse(zeros(n,1));       %Load vector
f(locf,:)=F;                %Assign loads
f(locsup,:)=[];             %Remove constrained DOFs
u=K\f;
B(locsup,:)=[];
dl = -(B'*u);       %Member elongation
N = E.*A./L.*dl;    %Member force
Sigma=N./A;         %Member stress

% Merging displacements with supports *********************
pr=1:n;
pr(locsup)=[];
U=sparse(zeros(n,1));
U(pr,1)=u;

p_deformed=[p(:,1)+U(1:2:n), p(:,2)+U(2:2:n)];

end


%% Volume as Objective Function *********************
function V = calculateVolume(L, A)

% Check that L and A have the same size
assert(length(L) == length(A), 'Length vectors L and A must have the same size');

% Calculate volume of each bar and sum
V = sum(L .* A);
end