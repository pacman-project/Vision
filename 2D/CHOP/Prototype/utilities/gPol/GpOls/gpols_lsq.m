function [mse,corcfsq,theta] = gpols_lsq(fss,X,Y);
%Calculates linear parameters and error values
%  [mse,corcfsq,theta,Yv] = gpols_lsq(fss,X,Y);
%    mse <- mean squared error
%    corcfsq <- square correlation coeffc.
%    theta <- identified linear parameters (last element is the bias)
%    fss -> cell array of function strings
%    X,Y -> regression matrix and output vector
%

% (c) Janos Madar, University of Veszprem, 2005

mse = Inf;
corcfsq = 0;
theta = zeros(length(fss)+1,1);

warning off

%Calculate model terms (eval functions)
XX = [];
for i = 1:length(fss),
  xxi = eval(fss{i});
  XX = [XX, xxi];
end
XX = [XX, ones(size(xxi))];

%Calculate parameters and estimated output
theta = XX \ Y;
YV = XX * theta;

warning on

%MSE:
mse = mean(sum((YV-Y).^2));
if isnan(mse), mse = Inf; end

%Correlation coefficient:
c = corrcoef(YV,Y);
c = c(1,2);
if isnan(c), c=0; end
corcfsq = max(0,min(c*c,1));
