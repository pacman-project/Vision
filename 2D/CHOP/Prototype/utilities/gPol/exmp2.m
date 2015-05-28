% Simple example for GP-OLS
%  Dynamical input-output model identification
% (Append the GpOls directory to the path)
%

clear all

%Simulation of a dynamic system and generates input/output data
t = [0:0.2:20]';
u = sin(t/2)-0.5;
u = u + randn(size(u))*0.1;
y = zeros(size(u));
y(1) = 0;
y(2) = 0;
for k = 3:length(t),
  dy = 0.7*u(k-1)*u(k-1) - 0.6*y(k-1) - 0.3*y(k-2) - 0.1;
  y(k) = y(k-1) + dy;
end
%Adds some simulated 'measurement noise' to the output
y = y + randn(size(y))*0.02;

%Select the maximum input and output order for identification
uorder = 2;
yorder = 2;

%Regressors and outputs for identification
tofs = max(uorder,yorder)+1;
Y = y(tofs:end) - y(tofs-1:end-1); %dy
X = [];
for  i=1:yorder,
  X = [X, y(tofs-i:end-i)];
end
for  i=1:uorder,
  X = [X, u(tofs-i:end-i)];
end

%GP equation symbols
symbols{1} = {'+','*'};
for i = 1:yorder,
  symbols{2}{i} = sprintf('y(k-%i)',i);
end
for j = 1:uorder,
  symbols{2}{i+j} = sprintf('u(k-%i)',j);
end

%Initial population
popusize = 40;
maxtreedepth = 5;
popu = gpols_init(popusize,maxtreedepth,symbols);

%first evaluation
opt = [0.8 0.7 0.3 2 2 0.2 25 0.01 1 0];
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

