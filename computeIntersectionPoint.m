function [intersectionCoords] = computeIntersectionPoint(coords1,angle1,coords2,angle2)
% This function compute the intersection point between two lines. 
% Inputs:
%   coords1 -- (x,y) of a point on the first line
%   angle1  -- angle of the first line with respect to x axis
%   coords2 -- (x,y) of a point on the second line
%   angle2  -- angle of the second line with respect to x axis
%
% Output:
%   intersectionCoords -- (x,y) of the intersection point between the two
%                         input lines


% Create set of dummy coordinates to populate the lines and find
% intersection. More coordinates is more accurate, and slower.
dummyCoordsX = linspace(-10,10,1000);

% Create line 1
line1Func = @(xInput) (xInput - coords1(1))*tand(angle1) + coords1(2);

% Create line 2
line2Func = @(xInput) (xInput - coords2(1))*tand(angle2) + coords2(2);

% Populate each line with the coordinates
line1X = dummyCoordsX;
line1Y = line1Func(dummyCoordsX);
line2X = dummyCoordsX;
line2Y = line2Func(dummyCoordsX);

% Find intersection point using polyxpoly();
[xInt,yInt] = polyxpoly(line1X,line1Y,line2X,line2Y);

intersectionCoords = [xInt,yInt];
end