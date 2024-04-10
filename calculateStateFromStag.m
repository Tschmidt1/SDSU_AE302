function [P,T,rho] = calculateStateFromStag(Ma,gamma,Po,To,rho_o)
% This function calculates Pressure (P), Temperature (T), and density (rho)
% using the isentropic flow relations. The inputs are Mach number (Ma),
% stagnation pressure (Po), stagnation temperature (To), specific heat ratio 
% (gamma), and stagnation density (rho_o);

%global gamma

% Pressure 
P = Po.*(1+(gamma-1)/2*Ma.^2).^(-gamma/(gamma-1));

% Temperature
T = To.*(1+(gamma-1)/2*Ma.^2).^-1;

% Density
rho = rho_o.*(1+(gamma-1)/2*Ma.^2).^(-1/(gamma-1));


end