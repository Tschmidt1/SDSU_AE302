function [y] = calculateNozzleContourFromEquation(AStar, x)
% This function computes a nozzle contour according to an equation. NOTE: this
% is different from the nozzle contour determined by the
% "ComputeMinimumNozzleContours.m". Ideally "ComputeMinimumNozzleContours.m"
% would be used to determine the (x,y) coordinates of the nozzle wall, then this
% function could be used to generate a spline. 
%
% In general, this function takes in x coordinates and computes the  y
% coordinates along a nozzle wall. 
%
% Inputs
%       AStar -- Nozzle throat area
%
%       x     -- x coordinates to compute nozzle wall y location at
%
% Outputs
%       y     -- y coordinates of the nozzle wall


% For simplicity, we use a simple y = x^(1/2) + yo equation to compute the
% nozzle wall.
nozzleWall = @(xIn,yIn) (xIn + 1.0).^(1/2) + yIn;

% Compute y coordinates
y = nozzleWall(x,AStar/2);

end