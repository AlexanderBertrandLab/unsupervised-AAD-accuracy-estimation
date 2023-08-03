function F = mMLE(mu,sigma,x)
% MMLE Maximum likelihood estimation function for the folded normal 
% distribution to solve.
%
%   Input parameters:
%       mu [DOUBLE]: mean of the underlying normal distribution
%       sigma [DOUBLE]: standard deviation of the underlying normal
%           distribution
%       x [DOUBLE]: vector of samples of the folded normal distribution
%
%   Output:
%       F [DOUBLE]: MLE-function evaluated at the given samples, with the
%           given parameters

% Author: Simon Geirnaert, KU Leuven, ESAT & Dept. of Neurosciences
% Correspondence: simon.geirnaert@esat.kuleuven.be

F = 0.5*sum(x-mu)-sum(x./(1+exp(2*mu.*x./sigma^2)));
