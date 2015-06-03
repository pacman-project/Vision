function [accuracyOverall] = Princeton_classification()

setenv OMP_NUM_THREADS 12
%addpath('libsvm-3.11/matlab');
addpath('categorization/SVM/SVM_chi2/libsvm_chi_openmp/libsvm_chi_ksirg/matlab');

threads = 12;

% this is structure to collect average error for every class
accuracies = 0;
accuracies = double(accuracies);
confMatr = zeros(20);

% assign category label to number


load('PrincetonTest.mat');
load('PrincetonTrain.mat');


% xx = load('xx_1Ch.mat');
% X = xx.X;
% Y = xx.Y;
% M = xx.M;

LabelsAll = [Ys_train; Ys_test];
numLabels = size(LabelsAll, 1);

CategoriesAll{1} = Ys_train{1};
cur = 1;

for i = 1:numLabels

    is_exist = false;
    
    for j = 1:size(CategoriesAll, 1);
        if strcmp(CategoriesAll{j}, LabelsAll{i});
            is_exist = true;
            break;
        end
    end
    
    if ~is_exist
        cur = cur+1;
        CategoriesAll{cur} = LabelsAll{i};
    end

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
    
    accuracies = accuracies + accuracy(1);
    % this is to compute confusion matrix   
    YY = [Ytest, predict_label];
    ll = length(Ytest);

    for i = 1:ll
        confMatr(YY(i,1), YY(i,2)) = confMatr(YY(i,1), YY(i,2)) + 1;
    end


accuracyOverall = accuracies 
 
% save('OneLeaveOut_confMatr1.mat', 'confMatr');
% save('OneLeaveOut_errors1.mat', 'errors1');
a = 2;


end



