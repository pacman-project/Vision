
% generate samples from Fischer von Mises distribution

n = 100; % number of points to be generated
mu = [0,0,1];
m = 3;
kappa = 10;

% V = rand(1,2);
% for i = 1:30 
%     W = 2 * rand(1) - 1;
%     X(i, :) = [sqrt(1-W^2)*V, W];
%     quiver3(0,0,0, X(i, 1), X(i,2), X(i, 3));
%     hold on
% end

X  = randvonMisesFisherm(m, n, kappa, mu);
X = X';

% for i = 1:n  
%     quiver3(0,0,0, X(i, 1), X(i,2), X(i, 3));
%     hold on
% end
% axis equal

xTrial = [0,1,0];
xTrial = xTrial/norm(xTrial);
vonMisesMultiplier = kappa/(2 * pi * (exp(kappa) - exp(-kappa)));

dist = 1 - (vonMisesMultiplier * exp(kappa*xTrial*mu')/kappa) / ((vonMisesMultiplier * exp(kappa*mu*mu')/kappa))


%% estimate parameters of fischer von mises given set of samples

% sumX = sum(X, 1);
% 
% muEst = sumX/norm(sumX)
% R = norm(sumX)/n;
% 
% kappaEst = R*(m - R^2)/(1-R^2) 



