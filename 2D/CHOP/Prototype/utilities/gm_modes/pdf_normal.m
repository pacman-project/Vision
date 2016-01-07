% [p,U,L,k] = pdf_normal(x,m,S[,U,L,k]) pdf of multivariate normal distribution
%
% In:
%   x: NxD list of D-dimensional row vectors.
%   m: Dx1 vector (mean).
%   S: DxD matrix (covariance).
% The following parameters are provided for efficiency (useful when
% one needs to apply pdf_normal to a lot of points for constant m, S).
%   U: DxD matrix with the eigenvectors of S columnwise (default: compute it).
%   L: Dx1 list of the eigenvalues of S (default: compute it).
%   k: normalisation constant (default: compute it).
% Out:
%   p: Nx1 list of values of the probability density p(x) at point x.
%   U: DxD matrix with the eigenvectors of S columnwise.
%   L: Dx1 list of the eigenvalues of S.
%   k: normalisation constant.

% Copyright (c) 1998 by Miguel A. Carreira-Perpinan

function [p,U,L,k] = pdf_normal(x,m,S,U,L,k)

% Argument defaults
if nargin==3
  [U,L] = eig(S);		% S = U*L*U'
  L = diag(L);
  k = 1/sqrt(prod(2*pi*L));
end;

p = k * exp(-sum(cdiv((U'*(x'-m(:,ones(1,length(x(:,1)))))).^2,L),1)/2)';

