function [Ma] = MachFromAOverAStar(supersonicStatus,AOverAStar,gamma)
% This function computes Mach number from A/A*
% Written by Professor Popov (SDSU)

% AOverAStar = 2; % A/A*
% 
% supersonicStatus = 0; % 1 for supersonic, 0 for subsonic
% 
% gamma = 1.2; % specific heat ratio


dM         = 1e-5; % small increment for numerical differentiation
iterations = 30; % how many times we want to iterate

if (supersonicStatus)
    Ma = 10; % if supersonic, start with a very high initial guess
else
    Ma = 0.1; % if subsonic, start with a very low initial guess
end

for ii = 1:iterations
    % evaluate A for the current M
    A1   = AOverAStarFromMach(Ma,gamma); 

    % evaluate the slope
    dAdM = (AOverAStarFromMach(Ma+dM,gamma) - AOverAStarFromMach(Ma-dM,gamma))./(2.*dM); 

    % take the Newton's method step
    Ma   = Ma + (AOverAStar - A1)./dAdM;                               
    
    %disp(['Mach number after iteration ',num2str(ii),' is: ',num2str(Ma)]); % output the result
end

end