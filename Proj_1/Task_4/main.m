%***********************************
%
%   Structural Optimization 2025
%   Olareanu Alexandru
%   Proj 1 task 4
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
stop = false;
i = 1;


while true
    Df = sensitivity(x,dx,fem,opts);

    %Stopping criteria
    if (norm(Df)<=1e-6) || (i>=30)
        % Stop optimization loop
        break
    end

    direction = -Df;
    % Inexact line search - Armijo's rule
    f_0 = Objective(x,fem,opts);
    f_0_prime = (Objective((x+dx*direction),fem,opts)-f_0)/dx; % FD aprox fro f in search direction
    alpha = 1;
    

    while true % Bracketing high side
        if Objective((x+alpha*eta*direction),fem,opts) > f_0 + alpha*rho*f_0_prime
            break
        else
            alpha = eta*alpha;
        end
    end


    alpha_hat = alpha;
    j = 1;
    while true % Armijo's rule
        alpha = mu^j * alpha_hat;
        if Objective((x+alpha*direction),fem,opts) <= f_0 + alpha*rho*f_0_prime
            x = x+alpha*direction;
            break
        else
            j = j+1;
        end
    end
    
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


