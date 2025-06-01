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
axis([-3 8 -1 12]);


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

% Visualize optimized structure
figure();
patch('vertices', p_deformed_opt, 'faces', t, 'facecolor', [.9, .9, .9]);
title('Deformed Mesh');
axis equal;
axis([-3 8 -1 12]);


% Print optimized cross-sectional areas and stresses
fprintf('\nOptimized cross-sectional areas and stresses:\n');
member_indices = (1:length(Sigma_opt))';
table(member_indices, A_optimal, Sigma_opt, N_opt, 'VariableNames', {'Member', 'Area', 'Stress', 'Force'})





% Find the node indices for visualization (convert from DOFs back to node numbers)
locsup_nodes = unique(ceil(locsup/2));
locf_node = unique(ceil(locf(1)/2));

% Plot initial structure
Plot_Structure(3, p, b, A_initial, 'Initial Structure - Volume:', sum(L.*A_initial), .1, 1e-3, locsup_nodes(1), locsup_nodes(2), locf_node);

% Plot optimized structure
Plot_Structure(4, p_deformed_opt, b, A_optimal, 'Optimized Structure - Volume:', fval, 5, 0.1, locsup_nodes(1), locsup_nodes(2), locf_node);





% Create new figure
newFig = figure;

% Define layout: 2 rows, 2 columns
rows = 2;
cols = 2;

% Loop through figures 1 to 4
for k = 1:4
    % Get handle to existing figure and its axes
    oldFig = figure(k);
    oldAx = findall(oldFig, 'type', 'axes');

    % Create subplot in new figure
    newAx = subplot(rows, cols, k, 'Parent', newFig);
    grid on;

    % Copy all axes objects to the new subplot
    for ax = flipud(oldAx')  % Flip to preserve visual order
        newCopiedAx = copyobj(ax, newFig);
        set(newCopiedAx, 'Position', get(newAx, 'Position'));
        delete(newAx);  % Delete placeholder subplot
    end
end












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