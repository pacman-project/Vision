% this is to perform recognition of all levels of the hierarchy

% do not forget: 'matlabpool open 4' for parallel computing

% dataSetNumber 
% for Aim@Shape datasetnumber = 1;
% for Washington datasetnumber = 2

function [] = LayersRecognition()

dataSetNumber = 2;
nClusters = 9;
n2Clusters = nClusters^2;

% define folders configuration
root = '/home/vvk201/LibHoP3D/';

addpath([root,'utility']);
addpath([root,'settings']);
addpath([root,'Optimization']);
addpath([root,'visualization']);
addpath([root,'statistics']);
addpath([root,'Recognition']);
addpath([root,'Learning']);



% define a path to the input files
if dataSetNumber == 1
    depthPathDefault = [root, 'settings\list_depth.mat'];
    depthPath = 'D:\3D\Input Data\Images for categorization\1T';    
elseif dataSetNumber == 2
    depthPath = '/home/vvk201/Wash-rgbd-dataset';
    depthPathDefault = '';
end

fileListPrecomputed = false;
is_subset = false; % whether we shall use all files for learning

% define the subset length
if dataSetNumber == 1
    subset_len = 400; % how much shall we use for training
    subsetPercent = 1.0; % not used
elseif dataSetNumber == 2
    subset_len = 1000;
    subsetPercent = 0.01; % what percent from each folder to use
end

% % downsampling (if necessary)
% if dataSetNumber == 1
%     is_downsampling = false;
%     dowsample_rate = 1;
% elseif dataSetNumber == 2    
%     is_downsampling = true;
%     dowsample_rate = 2.5;
% end

is_downsampling = false;
dowsample_rate = 1;

% define parameters sigma and sigmaKernelSize
if dataSetNumber == 1
    ss = load('settings/sigma.mat');
    sigma = ss.sigma; % this is a sigma for gaussian derivative 
    sigmaKernelSize = ss.sigmaKernelSize; % size of the kernel for gaussian
elseif dataSetNumber == 2
    sigma = 2.0;
    sigmaKernelSize = 9;
end 

displacement34 = 6;

% parameters for prelimiminary processing functions
isErrosion = false;
discSize = 1;
dxKernel = load('dxKernel.mat');
dxKernel = dxKernel.dxKernel;

% helpful stuff for an inhobition function
combs = load('settings/combs10.mat');
combs = combs.combs; % combinations for the line discretization function
largestLine = 10; % for decomposition of the line discretization function

% ------------ parameters for line discretization---------------------
% min [error + wOverlap * overlap - wCoverage*coverage ]
if dataSetNumber == 1
    wCoverage = 0.25;
    wOverlap = 0.4;
elseif dataSetNumber == 2    
    wCoverage = 3.25;
    wOverlap = 0.52;
end

% --------Define the layers for recognition--------------------------------
is_first_layer = false;
is_second_layer = false; 
is_third_layer = true;
is_4th_layer = false; 
is_5th_layer = false; 
is_6th_layer = false;

is_inhibition_first_layer = false;
is_inhibition_second_layer = false; 
is_inhibition_third_layer = false;
is_inhibition_4th_layer = false; 
is_inhibition_5th_layer = false; 
is_inhibition_6th_layer = false;

areLayersRequired = [is_first_layer, is_second_layer, is_third_layer, is_4th_layer, is_5th_layer, is_6th_layer];
isInhibitionRequired = [is_inhibition_first_layer, is_inhibition_second_layer, is_inhibition_third_layer, ...
                        is_inhibition_4th_layer, is_inhibition_5th_layer, is_inhibition_6th_layer];

% --------output folder-----------------------------------------------------


outRoot2 = '/home/vvk201/Wash-rgbd-dataset_layer2';
outRoot3 = '/home/vvk201/Wash-rgbd-dataset_layer3';
outRoot4 = '/home/vvk201/Wash-rgbd-dataset_layer4';
outRoot5 = '/home/vvk201/Wash-rgbd-dataset_layer5';
 
 
% --------input file names-----------------------------------------------------

dsN = num2str(dataSetNumber);
nCl = num2str(nClusters);
aL = '2'; % abstraction level

vocabulary1Layer = [root, 'statistics/statistics_1_', dsN, '_', nCl, '.mat'];
 
 % files for part selection results
parts3Layer = [root, 'statistics/partsSelectionResults_3_', dsN, '_', nCl, '_a', aL, '.mat'];
parts4Layer = [root, 'statistics/partsSelectionResults_4_', dsN, '_', nCl, '_a', aL, '.mat'];
parts5Layer = [root, 'statistics/partsSelectionResults_5_', dsN, '_', nCl, '_a', aL, '.mat'];
parts6Layer = [root, 'statistics/partsSelectionResults_6_', dsN, '_', nCl, '_a', aL, '.mat'];

% necessary files with statistics
% files for the Aggregated statistics
statistics1LayerAggregated = [root, 'statistics/statisticsAggregated_1_', dsN, '_', nCl, '.mat'];
statistics3LayerAggregated = [root, 'statistics/statisticsAggregated_3_', dsN, '_', nCl, '.mat'];
statistics4LayerAggregated = [root, 'statistics/statisticsAggregated_4_', dsN, '_', nCl, '.mat'];
statistics5LayerAggregated = [root, 'statistics/statisticsAggregated_5_', dsN, '_', nCl, '.mat'];
statistics6LayerAggregated = [root, 'statistics/statisticsAggregated_6_', dsN, '_', nCl, '.mat'];
 
 
% %--------------------------------------------------------------------------
% % creating a filelist here
% 
% if dataSetNumber == 1
%     [list_depth, lenF] = extractFileList(fileListPrecomputed, depthPath, depthPathDefault, is_subset, subset_len);
%     list_mask = [];
% elseif dataSetNumber == 2
%     [list_depth, list_mask, ~, lenF] = extractFileListWashington(fileListPrecomputed, depthPath, depthPathDefault, is_subset, subsetPercent);
% end
% 

%--------------------------------------------------------------------------

if is_first_layer || is_second_layer 
    
    % read the first layer vocabulary

    load(vocabulary1Layer); % should contain: cluster1Centres, cluster1Lengths, thresh
    
    computeCoverage(list_depth, list_mask, lenF, sigma, sigmaKernelSize, dxKernel, isErrosion, discSize, nClusters, ...
                vocabulary1Layer, areLayersRequired, outRoot2, combs, largestLine, wCoverage, wOverlap, isInhibitionRequired, ...
                is_downsampling, dowsample_rate, dataSetNumber, depthPath, cluster1Centres, cluster1Lengths, thresh);
end

% LIST_MASK - TO DELETE
            
if is_third_layer 
    
    elPath = outRoot2; % second layer elements
    is_subset = false;
    
    % [list_El, lenE] = extractFileList(fileListPrecomputed, elPath, elPath, is_subset, subset_len);
    [list_El, ~, ~, lenF] = extractFileListWashingtonForClassification(elPath, is_subset, 1.0);
    load(statistics3LayerAggregated); % 'X' ,'frequencies', 'curTS', 'triples'
    

    if dataSetNumber == 1
        load([root, 'statistics/layer3_WashBest.mat']);          % triples3Out
    elseif dataSetNumber == 2
        load(parts3Layer);   % partsOut, lenOut, abstractionLevel
    end

    lenOut = 600;
    partsOut = partsOut(1:lenOut, :);
    n3Clusters = lenOut;
    
    computeCoverage3_4(list_El, lenF, nClusters, n2Clusters, n3Clusters, X, partsOut, displacement34, abstractionLevel, ...
                areLayersRequired, outRoot3, combs, largestLine, wCoverage, wOverlap, isInhibitionRequired, is_downsampling, dowsample_rate, elPath);         
            
end


if is_4th_layer 
    elPath = outRoot3; % third layer elements
    
    if dataSetNumber == 1
        load('statistics/layer4WashBest.mat');                   % triples4Out
    elseif dataSetNumber == 2
        load(parts4Layer);   % partsOut
    end
    
    n4Clusters = size(partsOut, 1);
    
    computeCoverage_Vertical(list_el, lenF, nClusters, n2Clusters, n3Clusters, n4Clusters, triples3Out, triples4Out, displacement34, ...
            areLayersRequired, outRoot234, combs, largestLine, wCoverage, wOverlap, isInhibitionRequired, is_downsampling, dowsample_rate, elPath, ...
           atstractionFile3, abstractionFile4); 
    
end





