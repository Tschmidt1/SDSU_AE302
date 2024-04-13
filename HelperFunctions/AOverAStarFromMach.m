function output = AOverAStarFromMach(M,gamma)
% this is the area ratio for an isentropic nozzle, as a function of M and
% gamma
% Written by Professor Popov (SDSU)
output = (1./M).*( (  1 + ((gamma-1)./2).*M.^2  ) ./ ( (gamma+1)./2 )  ) .^ (  (gamma+1)./(2.*(gamma-1)) );