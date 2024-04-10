function [Po,To,rho_o] = calculateStagFromState(Ma,gamma,P,T,rho)
% This function calculates the stagnation properties of a flow at a given point
% using the isentropic flow relations. The inputs are the ratio of specific heats (gamma) 
% Mach number (Ma), pressure (P), temperature (T), and density (rho). The outputs
% are stagnation pressure (Po), stagnation temperature (To), and stagnation density (rho_o)

%global gamma

if exist('T','var')
    % Temperature
    To = T.*(1 + (gamma-1)/2*Ma.^2);
end

if exist('P','var')
    % Pressure 
    Po = P.*(1+(gamma-1)/2*Ma.^2).^(gamma/(gamma-1));
end

if exist('rho','var')
    % Density 
    rho_o = rho.*(1+(gamma-1)/2*Ma.^2).^(1/(gamma-1));
end

end