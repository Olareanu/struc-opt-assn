%***********************************
%
%   Structural Optimization 2025
%   Olareanu Alexandru
%   Proj 1 task Bonus
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



%% Task 3 dx estimation **************************


x_matrix=[-500 -500 ; 500 500; -500 500; 500 -500; 0 0]'; % x1 to x5

f_matrix = zeros(1,5);
for i=1:5
    f_matrix(i) = Objective(x_matrix(:,i)',fem,opts);
end

[max_error,i_max_error]=max(abs(f_matrix'-single(f_matrix')));
dx = double(sqrt(max_error / abs(Objective(x_matrix(:,i_max_error)',fem,opts))));


%% Task 4 Optimization

% initial point
x = [0,0];

% History arrays
xHistory = [];
FHistory = [];

% Initialize the console output with a header
disp('Iteration    Current x1    Current x2    Current f');

% hyperparameters
eta = 2;
rho = 0.2;
mu = 0.5;
% initialize
i = 1;

% Davidon's search, combination of slide 29, 28 and 35

% Please end this suffering 

G = eye(2);
Df = sensitivity(x, dx, fem, opts)';

while true
    % Stopping criteria
    if (norm(Df) <= 1e-6) || (i >= 30)
        % Stop optimization loop
        break
    end
    
    % Step 1
    direction = (-G * Df);
    f_0 = Objective(x, fem, opts);
    f_0_prime = (Objective((x + dx * direction'), fem, opts) - f_0) / dx;
    
    % Step 2
    alpha_1 = 0;
    f_1 = f_0;
    f_1_prime = f_0_prime;
    alpha_2 = 1;
    beta = 1;
    gamma = 0.9;
    x_old = x;
    f_2 = Objective(x + beta * direction', fem, opts);
    f_2_prime = (Objective((x + (beta + dx) * direction'), fem, opts) - f_2) / dx;
    
    while true
        if f_2_prime == 0
            alpha_star = beta;
            break
        end
        
        % Step 6
        if (f_2 > f_0) || (f_2_prime > 0)
            % Step 7
            A = [alpha_1^3      alpha_1^2       alpha_1 1;
                 3*alpha_1^2    2*alpha_1       1       0;
                 alpha_2^3      alpha_2^2       alpha_2 1;
                 3*alpha_2^2    2*alpha_2       1       0];
            fit_factors = [f_1;  f_1_prime;  f_2;    f_2_prime];
            cube = A\fit_factors;
            
            % Step 8
            alpha = (-cube(2) + sqrt(cube(2)^2 - 3*cube(1)*cube(3))) / (3*cube(1));
            f_alpha = Objective(x + alpha * direction', fem, opts);
            f_alpha_prime = (Objective((x + (alpha + dx) * direction'), fem, opts) - f_alpha) / dx;
            
            % Step 9
            if (f_alpha_prime >= gamma * f_0_prime) || (f_alpha_prime == 0)
                alpha_star = alpha;
                break
            end
            
            % Step 10
            if f_alpha_prime > 0
                alpha_1 = 0;
                alpha_2 = alpha;
                beta = alpha;
                f_2 = f_alpha;
                f_2_prime = f_alpha_prime;
            else
                alpha_1 = alpha;
                alpha_2 = beta;
                f_1 = f_alpha;
                f_1_prime = f_alpha_prime;
            end
        else
            beta = 2 * beta;
            alpha_2 = beta;
            f_2 = Objective(x + beta * direction', fem, opts);
            f_2_prime = (Objective((x + (beta + dx) * direction'), fem, opts) - f_2) / dx;
        end
    end
    
    x = x + alpha_star * direction';
    Df_old = Df;
    Df = sensitivity(x, dx, fem, opts)';
    gam = Df - Df_old;
    delta = (x - x_old)';
    G = G + (1 + (gam' * G * gam) / (delta' * gam)) * (delta * delta') / (delta' * gam) - ((delta * gam' * G + G * gam * delta') / (delta' * gam)); 
    
    % Store the current x and f in history arrays
    xHistory = [xHistory; x];
    FHistory = [FHistory; f_0];
    
    % Print the current state
    fprintf('%-11d\t%-12.4f\t%-12.4f\t%-12.4f\n', i, x(1), x(2), f_0);

    i = i + 1;
end

x_optimal = x;
f_optimal = Objective(x_optimal,fem,opts);

disp('x optimal:')
disp(x_optimal)
disp('f(x) optimal:')
disp(f_optimal)


% Plot History (as task 2)
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


