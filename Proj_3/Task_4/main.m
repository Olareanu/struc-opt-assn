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
warning('off', 'MATLAB:nearlySingularMatrix');
warning('off', 'MATLAB:singularMatrix');


[x_best, f_best, history] = simulatedAnnealingTruss();

function [x_best, f_best, history] = simulatedAnnealingTruss()
% Simulated Annealing for Truss Optimization
% Parameters from the image
k_max = 50;        % Number of iterations
m = 30;             % Number of moves per temperature
T0 = 20;           % Starting temperature
alpha = 0.85;      % Temperature reduction coefficient

% Measure execution time
tic;  % Start timer for entire algorithm

% Variable bounds
lb = [0.2, 0];   % Lower bounds for [h0, yf]
ub = [0.8, 2];   % Upper bounds for [h0, yf]

% Initialize
T = T0;
k = 1;

% Choose initial solution x
x(1) = 0.8;  % Initial height h0
x(2) = 0;  % Initial final height yf
x_best = x;

% Evaluate initial solution
[V_current, A_current] = optimizeTruss(x(1), x(2), false);
f_current = V_current;  % Assuming we want to minimize volume
f_best = f_current;

% History tracking - allocate for all possible evaluations
max_evaluations = k_max * m;
history.x = zeros(max_evaluations, 2);
history.f = zeros(k_max, 1);
history.T = zeros(k_max, 1);
history.x_best = zeros(k_max, 2);
history.f_best = zeros(k_max, 1);

eval_count = 0;  % Counter for all evaluations

fprintf('Starting Simulated Annealing...\n');
fprintf('Initial solution: h0=%.3f, yf=%.3f, f=%.6f\n', x(1), x(2), f_current);

% Main algorithm loop
while k <= k_max  % Stopping condition

    % Start timer for this iteration
    iter_tic = tic;

    % Inner loop: For i = 1 to m
    for i = 1:m

        % Generate candidate solution x_new in neighborhood of x
        x_new = generateNeighbor(x, lb, ub);

        % Evaluate new solution
        [V_new, A_new] = optimizeTruss(x_new(1), x_new(2), false);
        f_new = V_new;  % Objective function value

        % Store all evaluated points
        eval_count = eval_count + 1;
        history.x(eval_count, :) = x_new;

        % Calculate cost change
        delta_C = f_new - f_current;

        % Decision logic
        if delta_C <= 0
            % Accept improvement
            x = x_new;
            f_current = f_new;

            % Update best solution if better
            if f_new < f_best
                x_best = x_new;
                f_best = f_new;
            end

        else
            % Generate random number between 0 and 1
            randNr = rand();

            % Accept worse solution with probability exp(-delta_C/T)
            if exp(-delta_C / T) > randNr
                x = x_new;
                f_current = f_new;
            end
        end
    end

    % Store history
    history.x(k, :) = x;
    history.f(k) = f_current;
    history.T(k) = T;
    history.x_best(k, :) = x_best;
    history.f_best(k) = f_best;

    % Calculate iteration time
    iter_time = toc(iter_tic);

    % Print progress with iteration time
    fprintf('Iteration %d: T=%.3f, Current f=%.6f, Best f=%.6f, Time: %.2f sec\n', ...
        k, T, f_current, f_best, iter_time);

    % Temperature reduction: T^(k+1) = alpha * T^(k)
    T = alpha * T;

    % Use the best solution found at every subsequent step k as starting point
    x = x_best;
    f_current = f_best;

    % Increment iteration counter
    k = k + 1;
end

% Calculate total elapsed time
total_time = toc;

% Trim history.x to actual number of evaluations
history.x = history.x(1:eval_count, :);

fprintf('\nOptimization completed!\n');
fprintf('Best solution: h0=%.6f, yf=%.6f\n', x_best(1), x_best(2));
fprintf('Best objective value: %.6f\n', f_best);
fprintf('Total evaluations: %d\n', eval_count);
fprintf('Total time for Simulated Annealing: %.2f seconds\n', total_time);

% Plot results
plotResults(history, k_max);
optimizeTruss(x_best(1), x_best(2), true);

end




function x_new = generateNeighbor(x, lb, ub)
% Generate a neighbor solution within bounds

sigma = [0.05, 0.2]; % Different step sizes for h0 and yf
x_new = x + sigma .* randn(size(x));

% Ensure bounds are respected using clamping
x_new = max(lb, min(ub, x_new));

end




function plotResults(history, k_max)
% Plot optimization history

figure('Position', [100, 100, 1200, 800]);

% Plot objective function evolution
subplot(2, 2, 1);
plot(1:k_max, history.f, 'b-', 'LineWidth', 1);
hold on;
plot(1:k_max, history.f_best, 'r-', 'LineWidth', 2);
xlabel('Iteration');
ylabel('Objective Function Value');
title('Objective Function Evolution');
legend('Current', 'Best', 'Location', 'best');
grid on;

% Plot temperature evolution
subplot(2, 2, 2);
semilogy(1:k_max, history.T, 'g-', 'LineWidth', 2);
xlabel('Iteration');
ylabel('Temperature (log scale)');
title('Temperature Schedule');
grid on;

% Plot design variables evolution
subplot(2, 2, 3);
plot(1:k_max, history.x_best(:, 1), 'r-', 'LineWidth', 2);
hold on;
plot(1:k_max, history.x_best(:, 2), 'b-', 'LineWidth', 2);
xlabel('Iteration');
ylabel('Design Variable Value');
title('Design Variables Evolution');
legend('h0 (best)', 'yf (best)', 'Location', 'best');
grid on;

% Plot scatter plot of all evaluated design variables
subplot(2, 2, 4);
scatter(history.x(:, 1), history.x(:, 2), 50, 1:size(history.x, 1), 'filled');
colorbar;
xlabel('h0');
ylabel('yf');
title('Design Space Exploration (All Evaluations)');
grid on;
colormap(jet);
c = colorbar;
c.Label.String = 'Evaluation Order';
end

%% Truss Generation and Optimization

function [fval, A_optimal] = optimizeTruss(h0, yf, showPlots)

% Mesh generation (as per Exercise) *************
drectangle = @(p,x1,x2,y1,y2) -min(min(min(-y1+p(:,2),y2-p(:,2)),-x1+p(:,1)),x2-p(:,1));
fd = @(p)  drectangle(p,0,2,0,2);           %Domain
fh = @(p) ones(size(p,1),1);                % Uniform mesh
bounding_box=[0,0;2,2];                     %Domain bounding-box
p_keep=[0,0;0,2;2,0;2,2;0,0.5;0,1.5;2,yf];  %Points that are always kept in the generated mesh
[p,t,b,L] = distmeshSO( fd, fh, h0, bounding_box, p_keep );


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
options = optimoptions('fmincon', 'Display', 'off', 'Algorithm', 'sqp', 'MaxFunctionEvaluations', 1e5, ...
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

if showPlots == true
    figure(2);
    patch( 'vertices', p, 'faces', t, 'facecolor', [.9, .9, .9] )
    title('Initial Mesh');
    axis equal

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
    Plot_Structure(3, p, b, A_optimal, 'Optimized Structure (Undeformed) - Volume:', fval, 5, 0.1, locsup_nodes(1), locsup_nodes(2), locf_node);

    % Plot optimized structure for Load Case 1
    Plot_Structure(4, p_deformed_opt1, b, A_optimal, 'Optimized Structure (Deformed - Load Case 1) - Volume:', fval, 5, 0.1, locsup_nodes(1), locsup_nodes(2), locf_node);

    % Plot optimized structure for Load Case 2
    Plot_Structure(5, p_deformed_opt2, b, A_optimal, 'Optimized Structure (Deformed - Load Case 2) - Volume:', fval, 5, 0.1, locsup_nodes(1), locsup_nodes(2), locf_node);
end
end


%% Stress Constraint Function *********************
function [c, ceq] = stressConstraints(p, t, b, L, locf, locsup, A)
try

    % Calculate stresses for the current design
    [~, Sigma_1, ~, ~, Sigma_2, ~] = computeFe(p, t, b, L, locf, locsup, A);

    Sigma_y = 1;  % Yield stress

    % Inequality constraints: stress must be less than yield stress for both load cases
    c1 = abs(Sigma_1) - Sigma_y;
    c2 = abs(Sigma_2) - Sigma_y;
    c = [c1; c2]; % Combine constraints from both load cases

    % No equality constraints
    ceq = [];
catch ME
    warning('Error in constraint evaluation. Returning penalty values.');
    % Return large penalty values when constraint evaluation fails
    num_members = length(A);
    c = 1e6 * ones(2 * num_members, 1);  % Penalty for both load cases
    ceq = [];
end
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