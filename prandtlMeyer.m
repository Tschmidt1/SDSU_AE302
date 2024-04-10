function [pmf] = prandtlMeyer(Ma,gamma)
% This function calculates the Pranttl-Meyer function used in calculating flow
% pvariables across a Prandtl Meyer Expansion Fan. The inputs are Mach number
% (Ma) and specific heat ration (gamma)

pmf = sqrt((gamma+1)/(gamma-1))*atand(sqrt((gamma-1)/(gamma+1)*(Ma.^2-1))) - atand(sqrt(Ma.^2-1));
end