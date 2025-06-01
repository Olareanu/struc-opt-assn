%***********************************
%
%   Structural Optimization 2025
%   Olareanu Alexandru
%   Proj 3 task 1
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

%% Print mesh information
nNodes = size(p, 1);
nBars = size(b, 1);
barVolume = sum(L) * 0.1; % Total volume = sum of lengths * cross-sectional area

fprintf('Mesh information:\n');
fprintf('  Number of nodes: %d\n', nNodes);
fprintf('  Number of bars: %d\n', nBars);
fprintf('  Total bar volume: %.4f\n', barVolume);