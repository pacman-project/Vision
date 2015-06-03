% commonRoot = 'D:/'; 
% root = [commonRoot, 'LibHoP3D/'];
% addPaths(root);
% 
% addpath([root,'/yalmip']);
% addpath([root,'/yalmip/extras']);
% addpath([root,'/yalmip/demos']);
% addpath([root,'/yalmip/solvers']);
% addpath([root,'/yalmip/modules']);
% addpath([root,'/yalmip/modules/parametric']);
% addpath([root,'/yalmip/modules/moment']);
% addpath([root,'/yalmip/modules/global']);
% addpath([root,'/yalmip/modules/sos']);
% addpath([root,'/yalmip/operators']);

len = 100;

errs = rand(1,len)*2;
coverages = ones(1,len)*4;

A = diag(errs);
C = diag(coverages);

alpha = 1.2;

for i = 1:len
    if i ~= len
        A(i, i+1) = -0.25*errs(i);
        C(i, i+1) = -0.5* coverages(i);
    end
    if i ~= 1
        A(i, i-1) = -0.25*errs(i);
        C(i, i-1) = -0.5* coverages(i);
    end
end

x = binvar(len,1);

tic

obj = x'*(A-alpha*C)*x; 
ops = sdpsettings('solver','bnb','bnb.solver','fmincon');
sol = optimize([x(1:10) == 1,  integer(x)],obj,ops);

if sol.problem == 0
 % Extract and display value
 xx = value(x);
else
 display('Hmm, something went wrong!');
 sol.info
 yalmiperror(sol.problem)
end
xx'

toc;