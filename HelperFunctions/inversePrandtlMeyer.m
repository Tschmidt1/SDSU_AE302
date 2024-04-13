function [Ma] = inversePrandtlMeyer(nu)
% This function computes the inverse Prandtl-Meyer function: Ma(v). This is
% based on the paper: "Inverse Solutions of the Prandtl-Meyer Function" by
% Ozcan. This is an approximate solution as a general inversion of the PMF does
% not exist, but this approximate is found to be reasonably accurate in the
% paper referenced. 

global gamma

% maximum deflection angle for a gamma of 1.4; CHANGE TO MAKE GENERAL
nuMax = 130.45;     % degrees

% normalize nu
eta = nu/nuMax;

if (eta >= 0) && (eta <= 0.73)
    % Regime 1: 1 <= Ma < 7
    Ma = 1.03141 + 5.99196*eta - 15.6940*eta.^2 + 58.3798*eta.^3 ...
        - 87.6784*eta.^4 + 58.6275*eta.^5;
elseif (eta > .73) && (eta < 0.89)
    % Regime 2: 7 < Ma < 15
    Ma = 1429.53 - 7420.13*eta + 14512.2*eta.^2 - 12656*eta.^3 ...
        + 4169.33*eta.^4;
else
    disp('This function is not valid for Ma < 1 or Ma > 15!');
end


end

