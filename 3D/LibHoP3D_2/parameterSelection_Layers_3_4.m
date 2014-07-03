% this is to perform a grid search for different parameters of the third layer

% we have selected: bettaParams : 2.0
% gammaParam: 5.0 - DONE

% alphaParam = 0;
% betaParams = 1.6:0.4:6.0;    % 1.2  % 3
% gammaParams = 1.0:0.5:10.0;  % 5.0


% 4th layer
alphaParam = 0;
betaParams = 0.6 : 0.4: 1.8;    % 1.2  
gammaParams = [0.33, 0.66, 1];  % 5.0

% for the 4th layer we have selected
% we have selected: bettaParams : 1.2
% gammaParam: 0.66


root = 'D:\3D\LibHoP3D\';
addpath([root,'utility']);
addpath([root,'settings']);
addpath([root,'categorization']);
addpath([root,'Recognition']);
addpath([root,'Learning']);

n2Clusters = 81;

lenB = length(betaParams);
lenG = length(gammaParams);

accs = zeros(lenB, lenG);
num4Clusters = zeros(lenB, lenG);

for i = 1:lenB
    for j = 1:lenG
        [n3Clusters, n4Clusters] = learnHierarchy(alphaParam, betaParams(i), gammaParams(j));
        num4Clusters(i,j) = n4Clusters;
        if n4Clusters > 540
            continue;
        end
        LayersRecognition();
        histogramOfParts_best(n2Clusters, n3Clusters, n4Clusters);
        accuracyOverall = SVM_OneLeaveOut_classification();
        accs(i,j) = accuracyOverall;
    end
end


