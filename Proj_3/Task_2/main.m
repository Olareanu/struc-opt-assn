%***********************************
%
%   Structural Optimization 2025
%   Olareanu Alexandru
%   Proj 3 task 2
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
locf=locf*2-1; % x DOF loaded


%% FE ********************************************

% FE Parameters
F=.01;                      %Load
E=1;                        %Young's modulus
A=0.1*ones(size(b,1),1);    %Cross-section sizes

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

%% Visualizing deformed structure *********************
p_deformed=[p(:,1)+U(1:2:n), p(:,2)+U(2:2:n)];
figure();
patch( 'vertices', p_deformed, 'faces', t, 'facecolor', [.9, .9, .9] )
axis equal

%% Calculating reaction forces *********************
R=Kg*U;
fprintf('Reaction forces at at support nodes:\n')
table(locsup, R(locsup),'VariableNames',{'DOF', 'R'})

fprintf('Load forces at at load nodes:\n')
table(locf, R(locf),'VariableNames',{'DOF', 'F'})

%% Printing displacements at support and load nodes *********************
fprintf('\nDisplacements at support nodes:\n')
table(locsup, U(locsup),'VariableNames',{'DOF', 'Displacement'})

fprintf('\nDisplacements at load nodes:\n')
table(locf, U(locf),'VariableNames',{'DOF', 'Displacement'})


