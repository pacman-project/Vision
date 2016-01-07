% p = pdf_gm(x,mu,v,pm) pdf of isotropic Gaussian mixture
%
% Computes the value of p(x) for each point in x, where:
%   p(x) = \sum^M_{m=1}{p(m) p(x|m)}
% and p(x|m) is a Gaussian distribution of mean mu(m) and covariance v.
% Currently, all the mixture components must have the same, isotropic
% covariance (that is, v is a scalar independent of m).
%
% In:
%   x: NxD list of D-dimensional row vectors.
%   mu: MxD matrix (M centroids, rowwise).
%   v: positive scalar containing the variance (the same for each
%      component of the Gaussian mixture).
%   pm: Mx1 list (mixing proportions).
% Out:
%   p: Nx1 list of values of the probability density p(x) at point x.
%
% To do:
%
%   - Generalise to full-covariance mixture components, where v will
%     be a different covariance matrix for each component (currently it
%     assumes isotropic components, where v is a positive scalar, the
%     same for every component).
%
% See also pdf_normal, gm_moments, gm_modes_gq, gm_modes_fp.

% Copyright (c) 1999 by Miguel A. Carreira-Perpinan

function p = pdf_gm(x,mu,v,pm)

p = (exp(-sqdist(x,mu)/v/2)*pm)*(2*pi*v)^(-length(mu(1,:))/2);

