function [H,vH] = theirEntropy( hst )
%ENTROPY computes the entropy of a random variable (to base e)
%   [H,vH] = entropy( hst )
% INPUT:
%   hst -- histogram of the random variable values, e.g. the output of
%           hist(X,N) for some number of bins N.
% OUTPUT:
%   H -- real scalar -- estimated entropy in natural units
%   vH -- real scale -- standard error for the entropy estimate
%
% Implements eqn. (39) and (40) from:
%   "Estimating the errors on measured entropy and mutual information"
%   M.S.Roulston / Physica D 125 (1999) 285-294
%
% $Revision: 1.1 $
% By Shai Revzen, Berkeley 2007
hst = hst(:);
if any(hst<0) 
     error('entropy:negative','Negative entry in histogram data');
end
if all(hst==0)
    error('entropy:zero','Histogram is empty in all bins');
end
if any(hst~=round(hst))
    error('entropy:non-integer','Histogram must contain event counts -- i.e. integers');
end
% Normalize to probablity values
N = sum(hst);
hst = hst/N;
% Find support (non-empty bins)
spt = (hst>0);
% Compute ln(p) where meaningful
lh = zeros(size(hst));
lh(spt) = log(hst(spt));
% Observed entropy (eqn. (1))
Hobs = -sum(hst.*lh);
% Corrected entropy (eqn. (39))
Bs = length(hst)-sum(spt);
H = Hobs + (Bs-1)/(2*N);
% Standard error for the estimate (eqn. (40))
vH = [Hobs, Bs, N]; %sqrt( sum( (lh+Hobs).^2.*hst.*(1-hst))/ N );