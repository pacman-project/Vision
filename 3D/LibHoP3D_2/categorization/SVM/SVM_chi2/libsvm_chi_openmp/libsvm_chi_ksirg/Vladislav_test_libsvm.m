setenv OMP_NUM_THREADS 2

addpath('matlab');

tr_path='a1a.train';
tst_path='a1a.train';
 
[trYY trXX]=libsvmread(tr_path);
[tstYY tstXX]=libsvmread(tst_path);
trXXn=trXX;
%l1 - norm
trXXn=bsxfun(@rdivide,trXXn,sum(trXXn,2));
tic;
model = svmtrain(trYY, trXXn,'-c 4 -t 5 -z 2');
modelTime=toc;
 
tstXXn=tstXX;
tstXXn=bsxfun(@rdivide,tstXXn,sum(tstXXn,2)); %l1 - norm
 
tic
[pred, acc, dec_vals] = svmpredict(tstYY, tstXXn, model);
predTime = toc;
ss=sprintf('libsvm chi^2 acc=%0.5g modeltime=%g predtime=%g \n',acc(1),modelTime, predTime);
disp(ss);