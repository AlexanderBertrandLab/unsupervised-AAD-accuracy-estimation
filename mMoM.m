function F = mMoM(mu,sigma)
% MMOM Mean of the folded normal distribution.
%
%   Input parameters:
%       mu [DOUBLE]: mean of the underlying normal distribution
%       sigma [DOUBLE]: standard deviation of the underlying normal
%           distribution
%
%   Output:
%       F [DOUBLE]: mean of the folded normal distribution

% Author: Simon Geirnaert, KU Leuven, ESAT & Dept. of Neurosciences
% Correspondence: simon.geirnaert@esat.kuleuven.be

F = sigma*sqrt(2/pi)*exp(-mu^2/(2*sigma^2))+mu*erf(mu/sqrt(2*sigma^2));