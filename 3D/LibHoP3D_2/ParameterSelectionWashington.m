% this is a script to select the best set of parameters for washington
% dataset

scales = 1.5:0.5:5;
sigmas = 1:2:9;

lenScales = length(scales);
lenSigmas = length(sigmas);

% accs = zeros(lenScales, lenSigmas);
load('accuracies01.mat');
numChannels = 1;

for i = 7:8
    
    is_downsampling = false;
    dowsample_rate = scales(i);
    copyfile('D:/Input Data/Washington/Wash-rgbd-dataset_02', 'D:/Input Data/Washington/Wash-rgbd-dataset_02_scale');
    
    
    for j = 2:3
        
        curScale = scales(i);
        curSigma = sigmas(j);
        
        ratio = curSigma/curScale;
        
        if ratio > 5 || ratio < 1/4  % unfair ratio
            continue;
        end
            
        nNClusters = learnHierarchy(curSigma, is_downsampling, dowsample_rate); % learning and inference
        is_downsampling = false;
        
        histogramOfParts_best_Washington(nNClusters, numChannels);
        [accuracyOverall, confMatr] = SVM_OneLeaveOut_classification_Washington();

        accs(i,j) = accuracyOverall;
        
        save('accuracies01.mat', 'accs');
    end
    rmdir('D:/Input Data/Washington/Wash-rgbd-dataset_02_scale', 's');
    rmdir('D:/Input Data/Washington/Wash-rgbd-dataset_02_scale_layer2', 's');
end