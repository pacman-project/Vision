clear all;

setenv OMP_NUM_THREADS 4
%addpath('libsvm-3.11/matlab');
addpath('SVM_chi2/libsvm_chi_openmp/libsvm_chi_ksirg/matlab');

fold = 5;
threads = 4;

xx = load('xx.mat');
X = xx.X;
Y = xx.Y;
M = xx.M;

len = length(Y);

% convert to appropriate formats
nnn = isnan(X);
X(nnn) = 0;

% this is to split the set to training and test subsets
pTraining = 15;
pTest = 5;

MM = mod(M,20);
MM(MM == 0) = 20;

a = randperm(20);

indsTrain = [];
indsTest = [];

for i = 1:pTraining
    inds = find(MM == a(i));
    indsTrain = [indsTrain; inds];
end

for i = pTraining+1 : 20
    inds = find(MM == a(i));
    indsTest = [indsTest; inds];
end

Xtrain = X(indsTrain,:);
Xtest =  X(indsTest,:);
Ytrain = Y(indsTrain);
Ytest = Y(indsTest);

%[C, G, accuracy] = RadialCrossValidation(Xtrain, Ytrain, fold, threads);
%string = ['-t 2 -c ',num2str(C), ' -g ', num2str(G), ' -z 4'];  % Gaussian


[C, accuracy] = RadialCrossValidation_Chi(Xtrain, Ytrain, fold, threads);
string = ['-t 5 -c ',num2str(C), ' -z ', num2str(threads)];     % chi_squared


model = svmtrain(Ytrain, Xtrain, string);
[predict_label, accuracy, dec_values] = svmpredict(Ytest, Xtest, model);

% this is to compute confusion matrix

YY = [Ytest, predict_label];
load('confMatr.mat');

ll = length(Ytest);

for i = 1:ll
    confMatr(YY(i,1), YY(i,2)) = confMatr(YY(i,1), YY(i,2)) + 1;  
end

save('confMatr.mat', 'confMatr');





