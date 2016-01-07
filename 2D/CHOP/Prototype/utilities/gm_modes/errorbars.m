% [m, U, L] = errorbars(m,M,P[,f]) Confidence interval around a mode
%
% Given the mode of a distribution p(x), it computes error bars for the
% value of that mode in the form of a confidence interval (or a
% hyperrectangle, in several dimensions) containing a probability
% P. The unknown distribution is approximated around the mode m by a
% normal distribution having as mean m. The covariance matrix of this
% approximating normal is constructed in one of the following two ways:
%
% 1) So that the Hessian of the approximating normal is the same as
%    that of the unknown distribution. Use when the unknown
%    distribution is multimodal. This has two variants, depending on
%    whether the Hessian is from p(x) or from log(p(x)).
%
% 2) So that the covariance matrix of the approximating normal is the
%    same as that of the unknown distribution. Use when the unknown
%    distribution is unimodal.
%
% In any case, the appropriate matrix is contained in the formal
% parameter M. Which alternative to choose is decided by checking the
% sign of the matrix M (positive definite: covariance matrix, negative
% definite: Hessian). The parameter f tells, if M<0, whether it is
% from p(x) or from log(p(x)).
%
% This gives symmetric error bars, ignoring any skewness of the
% distribution around the mode.
%
% References:
%
%   Miguel A. Carreira-Perpinan (2000): "Mode-finding for mixtures of
%   Gaussian distributions", IEEE Trans. on Pattern Analysis and
%   Machine Intelligence 22(11): 1318-1323.
%
%   Miguel A. Carreira-Perpinan (1999): "Mode-finding for mixtures of
%   Gaussian distributions", Technical report CS-99-03, Dept. of
%   Computer Science, University of Sheffield, UK (revised Aug. 4, 2000).
%
% In:
%   m: 1xD vector containing the location of the mode.
%   M: DxD symmetric matrix containing either the Hessian of the
%      distribution at the mode (if M is negative definite) or the
%      covariance matrix of the distribution (if M is positive
%      definite).
%   P: real number in (0,1), the confidence level. Default: 0.9.
%   f: tells what function the Hessian corresponds to: 1 (to p(x)),
%      0 (to log(p(x))).
% Out:
%   m: 1xD vector containing the location of the mode (same as the input).
%   U: DxD orthogonal matrix containing the unit vectors of the
%      principal directions of the error bars.
%   L: 1xD vector containing the lengths of the error bars.
%
% The hyperrectangle centred in m, with sides aligned with the
% directions in U and with side lengths as given in L contains
% approximately a probability P (exactly for the approximating normal
% distribution).
%
% See also gm_modes_gq, gm_modes_fp.

% Copyright (c) 1999 by Miguel A. Carreira-Perpinan

function [m, U, L] = errorbars(m,M,P,f)

% ---------- Argument defaults ----------
if nargin==2 P=0.9; end;
if nargin<=3 f=0; end;

[U,L] = eig(M);			% U: orthonormal eigenvectors; L: eigenvalues
L = diag(L)';

D = length(L);
rho = sqrt(2)*erfinv(P^(1/D));

if max(L)<0			% M negative definite: Hessian
  L=-L;
  if f == 1
    L = prod(L/(2*pi))^(-1/(D+2))*L;
  end
  L = 2*rho*sqrt((ones(1,D)./L));
elseif min(L)>0			% M positive definite: covariance matrix
  L = 2*rho*sqrt(L);
else
  fprintf(stderr,'Matrix M is neither positive nor negative definite!\n');
end

