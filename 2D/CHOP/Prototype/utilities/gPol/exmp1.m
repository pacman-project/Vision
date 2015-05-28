% Simple example for GP-OLS
%  Static function identification
% (Append the GpOls directory to the path)
%

clear all

%Regression matrix
ndata = 100;
nvar = 3;
X = rand(ndata,nvar);

%Output vector (y = 10*x1*x2+5*x3)
Y = 10*X(:,1).*X(:,2) + 5*X(:,3);
Y = Y + randn(size(Y))*0.01; %some 'measurement' noise

%GP equation symbols
symbols{1} = {'+','*'};
symbols{2} = {'x1','x2','x3'};  %length(symbols{2}) = size(X,2) !

%Initial population
popusize = 40;
maxtreedepth = 5;
popu = gpols_init(popusize,maxtreedepth,symbols);

%first evaluation
opt = [0.8 0.7 0.3 2 1 0.2 30 0.05 0 0];
popu = gpols_evaluate(popu,[1:popusize],X,Y,[],opt(6:9));
%info
disp(gpols_result([],0));
disp(gpols_result(popu,1));
%GP loops
for c = 2:20,
  %iterate 
  popu = gpols_mainloop(popu,X,Y,[],opt);
  %info  
  disp(gpols_result(popu,1));
end

%Result
[s,tree] = gpols_result(popu,2);
disp(s);
