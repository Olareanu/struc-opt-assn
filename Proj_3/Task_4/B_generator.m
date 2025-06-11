function B = B_generator(p, b)
%***********************************
%   
%   Structural Optimization v2020
%   EDAC, ETHZ Zurich
%   Tino Stankovic
%   Yu Zhang
%   Exercise 2.2
%
%   Description: 
%      This function takes the following inputs: 
%          p:          nodes (represented by coordinates)
%          b:          bars  (each bar consists of starting node and ending node)     
%      It outputs two arrays as a representation of the ground structure:
%          B:          Full equilibrium matrix  
%***********************************
    %% Assign coordinates to nodes *************  
    B = zeros(size(p,1)*2,size(b,1)); % Initialize 
    %B = sym(zeros(size(p,1)*2,size(b,1))); % Initialize (comment Line 19 and uncomment this for symblic computation)            
    for idx_b = 1:size(b,1)       % Assign elements in B
       
        row1 = b(idx_b,1)*2;       % Row in B, starting node 
        x1 = (p(b(idx_b,2),1)-p(b(idx_b,1),1))/norm(p(b(idx_b,2),:)-p(b(idx_b,1),:));  % Calculate the trigonometric function
        y1 = (p(b(idx_b,2),2)-p(b(idx_b,1),2))/norm(p(b(idx_b,2),:)-p(b(idx_b,1),:));
        B(row1-1, idx_b) = x1;      
        B(row1, idx_b) = y1;

        row2 = b(idx_b,2)*2;       % Row in B, ending node 
        x2 = (p(b(idx_b,1),1)-p(b(idx_b,2),1))/norm(p(b(idx_b,1),:)-p(b(idx_b,2),:));  % Calculate the trigonometric function
        y2 = (p(b(idx_b,1),2)-p(b(idx_b,2),2))/norm(p(b(idx_b,1),:)-p(b(idx_b,2),:));
        B(row2-1, idx_b) = x2;
        B(row2, idx_b) = y2;
    end
end

