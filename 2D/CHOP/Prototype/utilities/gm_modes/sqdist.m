% sqd = sqdist(T[,Y,w]) Matrix of squared (weighted) distances between rows of T and Y
%
% If T and Y are matrices of row vectors of the same dimension (but
% possibly different number of rows), then D is a symmetric matrix
% containing the squared (weighted) Euclidean distances between every
% row vector in T with every row vector in Y.
%
% The square weighted Euclidean distance between vectors t and y is:
%   sqd = \sum_d { w_d (t_d - y_d)² }.
%
% In:
%   T: NxD matrix of N row D-dimensional vectors.
%   Y: KxD matrix of K row D-dimensional vectors (default equal to T).
%   w: 1xD vector of real numbers containing the weights (default: ones).
%
% Out:
%   sqd: NxK matrix of squared (weighted) Euclidean distances. D(n,k)
%   is the squared (weighted) Euclidean distance between row vectors
%   T(n,:) and Y(k,:).
%
% NOTE: this could be extended to compute the Mahalanobis distance
% given a positive definite matrix W: (t-y)'W(t-y).

% Copyright (c) 1999 by Miguel A. Carreira-Perpinan

function sqd = sqdist(T,Y,w)

% Argument defaults
if nargin==1 Y=T; end;

[N D] = size(T);
[K D] = size(Y);

if nargin>=3
  h = zeros(1,D); h(:) = sqrt(w);% We ensure that h is 1xD (even if w was Dx1).
  T = T.*h(ones(N,1),:);
  Y = Y.*h(ones(K,1),:);
end

% The intervector squared distance is computed as (t-y)² = t²+y²-2ty.
t = sum(T.^2,2);
y = sum(Y.^2,2);
sqd = t(:,ones(1,K)) + y(:,ones(N,1))' - 2*T*Y';

