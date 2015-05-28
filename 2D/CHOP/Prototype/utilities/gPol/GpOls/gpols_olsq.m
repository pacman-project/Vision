function [vfsdel,err] = gpols_olsq(fss,X,Y,limit);
% Calculates [err] values and decides which terms must be eliminated
%  [vfsdel,err] = gpols_olsq(fss,X,Y,limit)
%    vfsdel <- indexes of terms must be eliminated
%    err <- error reduction ratios
%    fss -> cell array of function strings
%    X,Y -> regression matrix and output vector
%    limit -> [err] treshold value
%

% (c) Janos Madar, University of Veszprem, 2005

vfsdel = [];
err = [];
if limit<=0,
  return;
end

%Calculate terms (eval functions)
warning off
XX = [];
for i = 1:length(fss),
  xxi = eval(fss{i});
  XX = [XX, xxi];
end
XX = XX - repmat(mean(XX),size(XX,1),1);
YY = Y - repmat(mean(Y),size(Y,1),1);
warning on

%OLS procedure
[err] = ols(XX,YY,0);
[dummy,ix] = sort(-err);

%Select:
if limit<1,
  %Select based on value
  for i = 2:length(ix),
    if err(ix(i)) < limit,
      vfsdel = [vfsdel, ix(i)];
    end
  end
else
  %Select based on first x pp.
  if limit<2, limit = 2; end
  for i = floor(limit):length(ix),
    vfsdel = [vfsdel, ix(i)];
  end
end 

%----------------------------------------------------------------------
function [err,theta] = ols(P,Y,q);
% OLS - Orthogonal Least Squares
%  [err,theta] = ols(P,Y,q);
%   P: regressor matrix (n*m) /m: number of regressors/
%   Y: output vector (n*1) /n: number of data points/
%   err: error reduction ratio vector (m*1)
%   theta: parameter vector (m*1)
%   q: q=0 -> [err] = ...; q~=0 -> [err,theta] = ...
%
%  Y = P*theta + e,
%   where Y: outputs, P: regressors, theta: parameters, e: error 

% (c) Janos Madar (Madár János), Univerisity of Veszprem, 2005

%Orthogonal decomposition ("economy size"):
[W,A] = qr(P,0);
%Error reduction ratios:
D = W'*W;
G = inv(D)*W'*Y;
yty = (Y'*Y);
err = zeros(size(P,2),1);
for i = 1:size(P,2),
  err(i) = G(i)*G(i)*W(:,i)'*W(:,i)/yty;
  %err(i) = (YY'*W(:,i))^2 / ((YY'*YY)*(W(:,i)'*W(:,i)));
end
%Parameter vector:
if q ~= 0,
  theta = inv(A)*G;
end
