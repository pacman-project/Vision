function [accuracyOverall, confMatr] = SVM_OneLeaveOut_classification_Washington()

tic
setenv OMP_NUM_THREADS 30
%addpath('libsvm-3.11/matlab');
addpath('categorization/SVM/SVM_chi2/libsvm_chi_openmp/libsvm_chi_ksirg/matlab');

threads = 30;

% this is structure to collect average error for every class
accuracies = [];
confMatr = zeros(51);

% X = [];
% Y = [];
% M = [];
% 
% for i = 1:2
%     str = ['xx', num2str(i),'.mat'];
%     X = [X; xx.X];
%     Y = [Y; xx.Y];
%     M = [M; xx.M]; % model is a number of the folder
% end

load('xx4.mat');
load('WashingtonTestSet.mat');  % testSet

len = length(Y);

% convert to appropriate formats
nnn = isnan(X);
X(nnn) = 0;

[numCat, numIter] = size(testSet);


for i = 1:numIter  % this is one which we leave out
    
    indsTest = [];
    indsTrain = [];
    
    for j = 1:numCat
        leaveOut = testSet(j,i);

        indsTestAdd =  find(M == leaveOut & Y == j);
        indsTrainAdd = find(M ~= leaveOut & Y == j);
        
        indsTest = [indsTest; indsTestAdd];
        indsTrain = [indsTrain; indsTrainAdd];
    end
    

    Xtrain = X(indsTrain,:);
    Xtest =  X(indsTest,:);
    Ytrain = Y(indsTrain);
    Ytest = Y(indsTest);

%   [C, G, accuracy] = RadialCrossValidation(Xtrain, Ytrain, fold);
%   string = ['-t 2 -c ',num2str(C), ' -g ', num2str(G)];

%     [C, accuracy] = RadialCrossValidation_Chi(Xtrain, Ytrain, fold, threads);
%     string = ['-t 5 -c ',num2str(C), ' -z ', num2str(threads)];    
%  no need in this cross-validation!!!
    
    string = ['-t 5 -c 5 -z ', num2str(threads)];
    
    model = svmtrain(Ytrain, Xtrain, string);
    [predict_label, accuracy, dec_values] = svmpredict(Ytest, Xtest, model);
    
    accuracies = [accuracies, accuracy(1)];
    % this is to compute confusion matrix   
    YY = [Ytest, predict_label];
    ll = length(Ytest);

    for i = 1:ll
        confMatr(YY(i,1), YY(i,2)) = confMatr(YY(i,1), YY(i,2)) + 1;
    end

end

accuracies = double(accuracies);
accuracies
accuracyOverall = sum(accuracies) / numIter;
 
% save('OneLeaveOut_confMatr1.mat', 'confMatr');
% save('OneLeaveOut_errors1.mat', 'errors1');

toc
end







