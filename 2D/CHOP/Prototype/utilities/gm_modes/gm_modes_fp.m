% [modes, pmodes, Hessians, codes, its] = 
%   gm_modes_fp(mu,v,pm[,tol,max_it,min_diff,max_eig])
% Modes of a Gaussian mixture by fixed-point iteration
%
% Finds all the modes of a Gaussian mixture, defined as:
%    p(x) = \sum^M_{m=1}{p(m) p(x|m)}
% where p(x|m) is a Gaussian distribution of mean mu(m) and covariance v.
% Currently, all the mixture components must have the same, isotropic
% covariance (that is, v is a scalar independent of m).
%
% gm_modes_fp searches for modes starting from every Gaussian centroid
% (contained rowwise in the matrix mu). This procedure should find all
% the modes. Repeated modes (due to the coalescence of several
% Gaussian components) are removed. The algorithm can be easily
% generalised to general mixtures of Gaussian distributions (not only
% the isotropic ones used here), but then there may be more modes than
% centroids.
%
% The mode search is done by a simple fixed-point iterative algorithm, 
% which needs no ad-hoc parameters to control convergence (e.g., step size).
%
% (Symmetric) error bars can be derived by approximating locally the
% mixture by a normal distribution with the same Hessian. The
% errorbars.m function does this.
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
%   mu: MxD matrix containing M D-dimensional centroids rowwise.
%   v: positive scalar containing the variance (the same for each
%      component of the Gaussian mixture).
%   pm: Mx1 list containing the mixing proportions of the mixture.
%   tol: minimum distance between successive points to keep iterating
%      (default 1e-4 standard deviations).
%   max_it: maximum number of iterations for the optimisation algorithm.
%      Default: 1000.
%   min_diff: minimum distance between modes to be considered the
%      same. Note that "tol" has to be small enough so that the
%      iterative procedure does not stop so early that the new mode
%      will not be considered different from the previously found
%      ones. Default: 1e-2 standard deviations.
%   max_eig: maximum algebraic value for an eigenvalue of the Hessian
%      to be considered definite negative. It is only used when
%      deciding whether to keep or throw away a fixed-point point, not
%      during the maximisation itself. Default: 0.
% Out:
%   modes: K x D matrix containing the K distinct modes rowwise. The
%      modes are returned sorted by descending probability (pmodes).
%   pmodes: K x 1 list containing the value of p(x) at each mode.
%   Hessians: D x D x K array. Hessians(:,:,k) is the Hessian at the 
%      kth mode.
%   codes: K x 1 list containing the exit code of the maximisation
%      algorithm (0: near-zero gradient found; 1: maximum number of
%      iterations reached).
%   its: K x 1 list containing the number of iterations employed in
%      the maximisation algorithm.
%
% To do:
%
%   - Generalise to full-covariance mixture components, where v will
%     be a different covariance matrix for each component (currently it
%     assumes isotropic components, where v is a positive scalar, the
%     same for every component).
%
% See also gm_modes_gq, gm_moments, errorbars.

% Copyright (c) 1999 by Miguel A. Carreira-Perpinan

function [modes, pmodes, Hessians, codes, its] = ...
    gm_modes_fp(mu,v,pm,tol,max_it,min_diff,max_eig)

% ---------- Argument defaults ----------
if nargin==3 tol=1e-4*sqrt(v); end;
if nargin<=4 max_it=1000; end;
if nargin<=5 min_diff=1e-2*sqrt(v); end;
if nargin<=6 max_eig=0; end;
% ---------- End of "argument defaults" ----------

[M,D] = size(mu);		% Number of components and dimensionality

% Discard all components whose mising proportion is smaller than a
% threshold (1% of the top value). This will speed up the search and
% probably won't miss any mode.
% This is heuristic. A worst case could be a lot of components with
% negligible probability but very close to each other, and an isolated 
% one with high probability. The cluster is discarded, as each
% individual component has negligible probability, but the combined
% probability for the cluster might be similar to that of the isolated 
% component. However, this is unlikely.
mu = mu(find(pm>max(pm)/100),:);
pm = pm(find(pm>max(pm)/100)); pm = pm/sum(pm);
M = length(mu(:,1));				% New number of components

% For mixtures of one component, the calculation is direct.
if M==1
  modes = mu;
  pmodes = (2*pi*v)^(-D/2);
  if nargout>=3
    Hessians = zeros(D,D,1); Hessians(:,:,1) = -pmodes*eye(D,D)/v;
  end
  if nargout>=4 codes = 0; end;
  if nargout>=5 its = 0; end;
  return;
end

modes = [];

for m=1:M	% Perform iterative procedure from every centroid in mu
  
  x = mu(m,:);					% Starting point: mth centroid
  x_old = x;
  code = -1; it = 0;
  % Compute the probability p(x). px_m and pxm are not multiplied
  % by the normalisation constant to save time and also to prevent
  % the numbers from getting too small [!]. I normalise in [1].
  while code < 0				% Fixed-point loop
    mux = mu - x(ones(M,1),:);			% mu_m - x for all m=1..M
    px_m = exp(-sum(mux.^2,2)/v/2);		% p(x|m) [!] for all m=1..M
    pxm = pm.*px_m;				% p(x,m|x) [!] forall m=1..M
    pxm = pxm*(2*pi*v)^(-D/2);			% Normalisation [1]
    px = sum(pxm);				% p(x)
    pm_x = pxm/px;				% p(m|x)
    x = pm_x'*mu;				% Fixed-point iteration
    
    if it>=max_it				% Max. no. iterations reached
      code = 1;
    elseif abs(x - x_old) < tol			% Tolerance achieved
      code = 0;
    else
      x_old = x;				% Keep current point
      it = it + 1;				% Continue iterating
    end
  end
  H = mux'*diag(pxm)*mux/v^2 - px*eye(D,D)/v;	% Hessian at x
  
  % ---------- Update mode list ----------
  % To avoid adding points with nonnegative Hessian (usually due to
  % areas of very low probability), we ensure that the maximum
  % eigenvalue (in algebraic value) is not positive as indicated by
  % max_eig.
  % Then we add (x,px,H) to the list if x is new, or if it was
  % already in the list but this one has higher probability, we
  % replace the old one(s). Note that it is necessary to replace not
  % just the old x in the list which is closest to the new one, but
  % we need to check every single old x in the list. For example, if
  % min_diff=0.1 and modes=[1.92 2.04] and x=2.01, the new list
  % should be modes=[2.01] and not [1.92 2.04].
  if max(eig(H)) < max_eig				% H<0 ?
    if isempty(modes)					% First mode found
      modes = x;
      pmodes = px;
      if nargout>=3 Hessians = zeros(D,D,1); Hessians(:,:,1) = H; end;
      if nargout>=4 codes = code; end;
      if nargout>=5 its = it; end;
    else
      temp1=sqrt(sum((modes-x(ones(length(modes(:,1)),1),:)).^2,2))>min_diff;
      same = find(temp1==0);		% These are all equivalent to x
      not_same = find(temp1==1);	% These are all different from x
      % Determine the old x with maximum probability
      if isempty(same)
	temp1=-1;
      else
	[temp1 temp2] = max(pmodes(same));
      end
      % Update modes, pmodes, Hessians, codes and its
      if temp1>=px					% Don't add the new x
	keep = [not_same;same(temp2)];
	modes = modes(keep,:);				% Update modes
	pmodes = pmodes(keep);				% Update pmodes
	if nargout>=3 Hessians=Hessians(:,:,keep); end;	% Update Hessians
	if nargout>=4 codes = codes(keep); end;		% Update codes
	if nargout>=5 its = its(keep); end;		% Update its
      else						% Add the new x
	modes = [modes(not_same,:); x];			% Update modes
	pmodes = [pmodes(not_same); px];		% Update pmodes
	if nargout>=3					% Update Hessians
	  temp1 = zeros(D,D,length(pmodes));
	  if ~isempty(not_same)
	    temp1(:,:,1:length(pmodes)-1) = Hessians(:,:,not_same);
	  end
	  temp1(:,:,length(pmodes)) = H;
	  Hessians = temp1;
	end
	if nargout>=4 codes=[codes(not_same);code]; end;% Update codes
	if nargout>=5 its = [its(not_same); it]; end;	% Update its
      end
    end
  end
  % ---------- End of "update mode list" ----------

end

% ---------- Sort modes by descending probability ----------
[temp1 temp2] = sort(-pmodes);
modes = modes(temp2,:);
pmodes = pmodes(temp2);
if nargout>=3 Hessians = Hessians(:,:,temp2); end;
if nargout>=4 codes = codes(temp2); end;
if nargout>=5 its = its(temp2); end;
% ---------- End of "sort modes by descending probability" ----------

