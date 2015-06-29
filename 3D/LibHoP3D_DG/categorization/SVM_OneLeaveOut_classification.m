function [accuracyOverall] = SVM_OneLeaveOut_classification()

setenv OMP_NUM_THREADS 12
%addpath('libsvm-3.11/matlab');
addpath('categorization/SVM/SVM_chi2/libsvm_chi_openmp/libsvm_chi_ksirg/matlab');

threads = 12;

% this is structure to collect average error for every class
accuracies = 0;
accuracies = double(accuracies);
confMatr = zeros(20);

xx = load('xx.mat');
X = xx.X;
Y = xx.Y;
M = xx.M;

len = length(Y);

% convert to appropriate formats
nnn = isnan(X);
X(nnn) = 0;

% this is to split the set to training and test subsets
pTraining = 19;
pTest = 1;
MM = mod(M,20);
MM(MM == 0) = 20;

% compute a matrix A for cross-bin similarities
parts = zeros(81,2);
for i = 1:81
    [parts(i,1), parts(i,2)] = compute2derivatives(i, 9);
end

A = Isodata_distances(parts, parts, 81, 81, false, false);

linkQC = @QC;


for iii = 1:5  % this is one which we leave out
    leaveOut = iii;

    indsTrain = [];
    indsTest = [];
    
    a1 = 1:leaveOut - 1;
    a2 = leaveOut+1:20;
    a3 = [leaveOut];
    a = [a1,a2,a3];
    

    for i = 1:19
        inds = find(MM == a(i));
        indsTrain = [indsTrain; inds];
    end

    for i = 20 : 20
        inds = find(MM == a(i));
        indsTest = [indsTest; inds];
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
    
    % here we compute a kernel function
    
    
    
    % string = ['-t 6 -c 5 -z ', num2str(threads)];
%     model = svmtrain(Ytrain, Xtrain, string);
%     [predict_label, accuracy, dec_values] = svmpredict(Ytest, Xtest, model);
    
    K =  [ (1:numTrain)' , rbfKernel(trainData,trainData) ];
    KK = [ (1:numTest)'  , rbfKernel(testData,trainData)  ];

    string = ['-t 4 -z ', num2str(threads)]; 
    model = svmtrain(Ytrain, K, string);
    [predict_label, accuracy, dec_values] = svmpredict(Ytest, KK, model);
    
    accuracies = accuracies + accuracy(1);
    % this is to compute confusion matrix   
    YY = [Ytest, predict_label];
    ll = length(Ytest);

    for i = 1:ll
        confMatr(YY(i,1), YY(i,2)) = confMatr(YY(i,1), YY(i,2)) + 1;
    end

end

accuracyOverall = accuracies / 20
 
% save('OneLeaveOut_confMatr1.mat', 'confMatr');
% save('OneLeaveOut_errors1.mat', 'errors1');
a = 2;


% Embedded function for cross-bin similarities
function dist = QC(P,Q,A,m)
    Z= (P+Q)*A;
    % 1 can be any number as Z_i==0 iff D_i=0
    Z(Z==0)= 1;
    Z= Z.^m;
    D= (P-Q)./Z;
    % max is redundant if A is positive-semidefinite
    dist= sqrt( max(D*A*D', 0));
end

end


