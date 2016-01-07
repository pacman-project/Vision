% [Mean, C, pMean] = gm_moments(mu,v,pm) Moments of a Gaussian mixture
%
% Computes the mean and the covariance matrix of a Gaussian mixture,
% defined as:
%    p(x) = \sum^M_{m=1}{p(m) p(x|m)}
% where p(x|m) is a Gaussian distribution of mean mu(m) and covariance v.
% Currently, all the mixture components must have the same, isotropic
% covariance (that is, v is a scalar independent of m).
%
% Actually this function is valid for any mixture, not just Gaussian.
%
% In:
%   mu: MxD matrix containing M D-dimensional centroids rowwise.
%   v: positive scalar containing the variance (the same for each
%      component of the Gaussian mixture).
%   pm: Mx1 list containing the mixing proportions of the mixture.
% Out:
%   Mean: 1xD vector containing the mean.
%   C: DxD symmetric positive definite matrix containing the covariance.
%   pMean: real number containing the value of p(x) at the mean.
%
% To do:
%
%   - Generalise to full-covariance mixture components, where v will
%     be a different covariance matrix for each component (currently it
%     assumes isotropic components, where v is a positive scalar, the
%     same for every component).
%
% See also gm_modes_gq, gm_modes_fp, errorbars.

% Copyright (c) 1999 by Miguel A. Carreira-Perpinan

function [Mean, C, pMean] = gm_moments(mu,v,pm)

[M,D] = size(mu);		% Number of components and dimensionality

Mean = pm'*mu;
C = v*eye(D,D) - Mean'*Mean;
for m=1:M
  C = C + pm(m)*mu(m,:)'*mu(m,:);
end

pMean = pdf_gm(Mean,mu,v,pm);

