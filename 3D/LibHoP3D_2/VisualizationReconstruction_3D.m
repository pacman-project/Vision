% this is for 3D visualization of the objects represented in terms of
% vocabulary

dataSetNumber = 2;
nClusters = 7;

% define folders configuration
commonRoot = 'D:/'; 
root = [commonRoot, 'LibHoP3D/'];
addPaths(root);

dsN = num2str(dataSetNumber);
nCl = num2str(nClusters);
vocabulary1Layer = [root, 'statistics/statistics_1_', dsN, '_', nCl, '.mat'];

load(vocabulary1Layer);

% define a path to the input files

depthPathDefault = '';
if dataSetNumber == 1
    depthPathDefault = [root, 'settings/list_depth.mat'];
    depthPath = [commonRoot, 'Input Data/AimShape/4T_600'];    
elseif dataSetNumber == 2
    depthPath = [commonRoot,'Input Data/Washington/Wash-rgbd-dataset_001'];    
elseif dataSetNumber == 3
    depthPath = [commonRoot, 'Input Data/VladislavSTD/Vladislav_STD/depth'];     
end

% this are lists of files with detections:
outRoot2 = [depthPath, '_layer2_inhibition'];
outRoot3 = [depthPath, '_layer3_inhibition'];
outRoot4 = [depthPath, '_layer4_inhibition'];
outRoot5 = [depthPath, '_layer5_inhibition'];
outRoot6 = [depthPath, '_layer6_inhibition'];
outRoot7 = [depthPath, '_layer7_inhibition'];
outRoot8 = [depthPath, '_layer8_inhibition'];

fileListPrecomputed = false;
is_subset = false; % whether we shall use all files for learning

% define the subset length
if dataSetNumber == 1
    subset_len = 400; % how much shall we use for training
    subsetPercent = 1.0; % not used
elseif dataSetNumber == 2
    subset_len = 1000;
    subsetPercent = 0.01; % what percent from each folder to use
elseif dataSetNumber == 3
    subset_len = 1000;
    subsetPercent = 0.01; % what percent from each folder to use
end

if dataSetNumber == 1 || dataSetNumber == 3
    [list_depth, lenF] = extractFileList(fileListPrecomputed, depthPath, depthPathDefault, is_subset, subset_len);
    [list_El, lenF] = extractFileList(fileListPrecomputed, elPath, '', is_subset, subset_len);
    list_mask = [];
elseif dataSetNumber == 2
    [list_depth, list_mask, ~, lenF] = extractFileListWashington(fileListPrecomputed, depthPath, depthPathDefault, is_subset, subsetPercent);
    [list_els, ~, ~, lenF] = extractFileListWashingtonForClassification(outRoot2, is_subset, subsetPercent);
end

[dxKernel, combs, largestLine, sigma, sigmaKernelSize, isErrosion, discRadius, is_guided, r_guided, eps, ...
    is_mask_extended, maxExtThresh1, maxExtThresh2] = loadFilteringParameters(dataSetNumber);
isY = true;
isX = true;
isTrim = false;

is_first_layer = false;
is_second_layer = true;
is_third_layer = false;
is_4th_layer = false; 
is_5th_layer = false; 
is_6th_layer = false;


depthStep = thresh/4;
offsetD = 20;    

for i = 1:1 % lenF
    I = imread(list_depth{i});
    marks = imread(list_els{i});
    
    [I, ~, ~, ~, r, c, is_successfull] = preliminaryProcessing(I, [], isErrosion, discRadius, isX, isY, ...
                            isTrim, dxKernel, sigmaKernelSize, sigma, is_guided, r_guided, eps, is_mask_extended, maxExtThresh1, maxExtThresh2);
    
%     imtool(I);
%     imtool(marks);
    
    maxRC = max([r,c]);
    [rows, cols] = find(marks>0);
    numEl = length(rows);
    depths = zeros(1,numEl);
    
    for j = 1:numEl
        % extract depths
        depths(j) = I(rows(j), cols(j));
    end
    
    
    minD = min(depths);
    maxD = max(depths);
    depthRange = ceil((maxD - minD)/depthStep) + 2 * offsetD;
    f = figure;
    
    fieldSize = [maxRC, maxRC, depthRange];
    
    for j = 1:numEl
        elements = marks(rows(j), cols(j));
        curDepth = floor((I(rows(j), cols(j)) - minD ) / depthStep) + offsetD;
        positions = [cols(j),rows(j),curDepth];
        
        
        % surfaceVisualizer(fieldSize, positions, elements, nClusters, cluster1Centres, depthStep, faceColor);
        [out] = surfaceVisualizer(fieldSize, positions, elements, nClusters, cluster1Centres, depthStep);
        hold on
    end

 
%     surfaceVisualizer(fieldSize, positions, 1, nClusters, cluster1Centres, depthStep);

    
    hold off
end



