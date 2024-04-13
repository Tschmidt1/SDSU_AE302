function [M2, R2, P2, T2] = calculateNormalShock_Downstream(M1, R1, P1, T1, gamma, R)
% This function calculates the flow properties downstream of a normal shock.
% The derivation of these equations can be found in John. D. Anderson's "Modern
% Compressible Flow," Chapter 3.
%
% Inputs:
%       M1    -- Upstream Mach number
%
%       R1    -- Upstream density
%
%       P1    -- Upstream pressure
%
%       T1    -- Upstream temperature
%
%       gamma -- ratio of the specific heats (Cp/Cv)
%
%       R     -- gas constant
%
% Outputs:
%       M2    -- Downstream Mach number
%
%       R2    -- Downstream density
%
%       P2    -- Downstream pressure
%
%       T2    -- Downstream temperature

% Create variables for common values
gammaMinus1    = gamma - 1;
gammaPlus1     = gamma + 1;
M1Squared      = M1.^2;
gammaM1Squared = gamma .* M1Squared;

M2 = sqrt(  ( 1 + gammaMinus1./2 .*M1Squared )./( gammaM1Squared - gammaMinus1./2 )  );
R2 = R1.*( (gammaPlus1.*M1Squared)./(2 + gammaMinus1.*M1Squared) );
P2 = P1.*( 1 + (2.*(gamma)./gammaPlus1).*(M1Squared - 1) );
T2 = P2./(R.*R2);

end