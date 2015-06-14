% this is to try a pre-defined kernel
function PreDefinedKernelTrial()

    setenv OMP_NUM_THREADS 30
    
    threads = 30;
    %addpath('libsvm-3.11/matlab');
    addpath('categorization/SVM/SVM_chi2/libsvm_chi_openmp/libsvm_chi_ksirg/matlab');

    is_GPU_USED = true;

    aa = load('./xx4_2.mat');
    X = aa.X;
    Y = aa.Y;
    M = aa.M;
    clear('aa');
    
    nNClusters{1} = 7;
    nNClusters{2} = 49;

    XX = (1:49)';
    lenCombs = 49;
    layerID = 2;
    
    displ3 = 6;
    displ5 = 18;
    displ7 = 52;
    fieldSize = [17,17,71];
    fieldCenter = ceil(fieldSize/2);
    
    fileForVisualizationPrevLayer{layerID-1} = '';
    downsamplingScheme = 3;
    clusterCurDepths = [];
   
    aa = load('statistics\statistics_1_2_7.mat');  % cluster1Centres
    cluster1Centres = aa.cluster1Centres;
    clear('aa');
    
    aa = load('WashingtonTestSet.mat');  % testSet
    testSet = aa.testSet;
    clear('aa');
   
    [X_first, nPrevClusters, emptyIndicator] = convertToSurfaceDescriptor(XX, lenCombs, layerID, nNClusters{1}, nNClusters{2}, fileForVisualizationPrevLayer{layerID-1},  ...
                       displ3, displ5, displ7, fieldCenter, cluster1Centres, downsamplingScheme, clusterCurDepths);  % this is a first layer descriptor

    AA = Integral_distances(X_first, X_first, nNClusters{2}, nNClusters{2}, false, false);
    
    AA = exp(-AA/0.5);
    A = zeros(15*49+4);
    
    for i = 1:15 % for each bin
        curOffset = (i-1)*49 + 1; 
        A(curOffset:curOffset+49-1, curOffset:curOffset+49-1) = AA;
    end
    
    m = 0.5;
    accuracies = [];

    [numCat, numIter] = size(testSet);
    for i = 1:1%numIter  % this is one which we leave out

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
        
        % compute kernels (distances between histograms): -----------------
        
        if is_GPU_USED
            Xtrain = gpuArray(Xtrain);
            Xtest = gpuArray(Xtest);
            A = gpuArray(A);
        end
        
        
%         K = ComputeQuadraticChiKernel(Xtrain, Xtrain, A, m, is_GPU_USED);
%         save('Temp/K.mat', 'K', '-v7.3');
%         KK = ComputeQuadraticChiKernel(Xtest,  Xtrain, A, m, is_GPU_USED);
%         save('Temp/KK.mat', 'KK', '-v7.3');

        K = ComputeChi_Squared_Kernels(Xtrain, Xtrain, is_GPU_USED);
        KK = ComputeChi_Squared_Kernels(Xtest, Xtrain, is_GPU_USED);


        a = load('Temp/KK.mat');
        KK = a.KK;
        clear('a');
        a = load('Temp/K.mat');
        K = a.K;
        clear('a');
        
        %------------------------------------------------------------------
        
        K =  [ (1:length(Xtrain))', K] ;  %  rbfKernel(Xtrain,Xtrain)
        KK = [ (1:length(Xtest))', KK];  %  rbfKernel(Xtest,Xtrain) 

        string = ['-t 4 -z ', num2str(threads)];
        
        model = svmtrain(Ytrain, K, string);
        [predict_label, accuracy, dec_values] = svmpredict(Ytest, KK, model);
        
%        string = ['-t 5 -c 5 -z ', num2str(threads)];
%         model = svmtrain(Ytrain, Xtrain, string);
%         [predict_label, accuracy, dec_values] = svmpredict(Ytest, Xtest, model);

        accuracies = [accuracies, accuracy(1)];
        % this is to compute confusion matrix   
        YY = [Ytest, predict_label];
        ll = length(Ytest);

%         for ii = 1:ll
%             confMatr(YY(ii,1), YY(ii,2)) = confMatr(YY(ii,1), YY(ii,2)) + 1;
%         end

    end

    accuracies = double(accuracies);
    accuracies
    accuracyOverall = sum(accuracies) / numIter;
    





% the following is an example code!

%# read dataset
%[dataClass, data] = libsvmread('C:\Projects\Vladislav\LibHoP3D\categorization\SVM\libsvm-2.9-dense_chi_square_mat/heart_scale');
% %# split into train/test datasets
% trainData = data(1:150,:);
% testData = data(151:270,:);
% trainClass = dataClass(1:150,:);
% testClass = dataClass(151:270,:);
% numTrain = size(trainData,1);
% numTest = size(testData,1);
% 
% %# radial basis function: exp(-gamma*|u-v|^2)
% sigma = 2e-3;
% rbfKernel = @(X,Y) exp(-sigma .* pdist2(X,Y,'euclidean').^2);
% 
% %# compute kernel matrices between every pairs of (train,train) and
% %# (test,train) instances and include sample serial number as first column
% K =  [ (1:numTrain)' , rbfKernel(trainData,trainData) ];
% KK = [ (1:numTest)'  , rbfKernel(testData,trainData)  ];
% 
% %# train and test
% model = svmtrain(trainClass, K, '-t 4');
% [predClass, acc, decVals] = svmpredict(testClass, KK, model);
% 
% %# confusion matrix
% C = confusionmat(testClass,predClass);

end
