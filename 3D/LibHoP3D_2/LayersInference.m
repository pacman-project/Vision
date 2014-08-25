% this is to perform inference of all levels of the hierarchy

% do not forget: 'matlabpool open 4' for parallel computing

% dataSetNumber 
% for Aim@Shape datasetnumber = 1;
% for Washington datasetnumber = 2

function [] = LayersInference()

dataSetNumber = 2;
nClusters = 7;
n2Clusters = nClusters^2;

% % here we initialize the perallel computing
% g = gpuDevice(1);
% reset(g);

% matlabpool open 10    % for parallel computin

% define folders configuration
commonRoot = 'D:/'; 
root = [commonRoot, 'LibHoP3D/'];
addPaths(root);


% define a path to the input files

depthPathDefault = 'listDepthDefault.mat';
if dataSetNumber == 1
    depthPathDefault = [root, 'settings/list_depth.mat'];
    depthPath = [commonRoot, 'Input Data/AimShape/4T_600'];    
elseif dataSetNumber == 2
    depthPathDefault = 'settings/listDepthDefault.mat';
    depthPath = [commonRoot,'Input Data/Washington/Washington3Categories_008'];      % '/home/vvk201/Wash-rgbd-dataset_0003T' 
elseif dataSetNumber == 3
    depthPath = [commonRoot, 'Input Data/VladislavSTD/Vladislav_STD/depth'];     
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
elseif dataSetNumber == 3
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

% define all filtering parameters 

[dxKernel, combs, largestLine, sigma, sigmaKernelSize, isErrosion, discRadius, is_guided, r_guided, eps, ...
    is_mask_extended, maxExtThresh1, maxExtThresh2] = loadFilteringParameters(dataSetNumber);

% ------------ parameters for line discretization---------------------
% min [error + wOverlap * overlap - wCoverage*coverage ]
if dataSetNumber == 1 || dataSetNumber == 3
    wCoverage = 0.25;
    wOverlap = 0.4;
elseif dataSetNumber == 2    
    wCoverage = 0.05;
    wOverlap = 0.05;
end


[ meargeThresh3, meargeThresh4,meargeThresh5, meargeThresh6, meargeThresh7, meargeThresh8 ] = defineMeargeThreshes(1);

% --------Define the layers for inference--------------------------------
is_first_layer = false;
is_second_layer = false;
is_third_layer = false;
is_4th_layer = false; 
is_5th_layer = false; 
is_6th_layer = true;
is_7th_layer = false; 
is_8th_layer = false;

is_inhibition_first_layer = false;
is_inhibition_second_layer = false; 
is_inhibition_third_layer = false;
is_inhibition_4th_layer = false; 
is_inhibition_5th_layer = false; 
is_inhibition_6th_layer = false;
is_inhibition_7th_layer = false; 
is_inhibition_8th_layer = false;


is_reconstructionError_first_layer = true;
is_reconstructionError_second_layer = false; 
is_reconstructionError_third_layer = false;
is_reconstructionError_4th_layer = false; 
is_reconstructionError_5th_layer = false; 
is_reconstructionError_6th_layer = false;
is_reconstructionError_7th_layer = false; 
is_reconstructionError_8th_layer = false;

areLayersRequired = [is_first_layer, is_second_layer, is_third_layer, is_4th_layer, is_5th_layer, is_6th_layer, is_7th_layer, is_8th_layer];
isInhibitionRequired = [is_inhibition_first_layer, is_inhibition_second_layer, is_inhibition_third_layer, ...
                        is_inhibition_4th_layer, is_inhibition_5th_layer, is_inhibition_6th_layer, is_inhibition_7th_layer, is_inhibition_8th_layer];
isReconstructionErrorRequired = [is_reconstructionError_first_layer, is_reconstructionError_second_layer, is_reconstructionError_third_layer, ...
    is_reconstructionError_4th_layer, is_reconstructionError_5th_layer, is_reconstructionError_6th_layer, is_reconstructionError_7th_layer, is_reconstructionError_8th_layer];

% --------output folder-----------------------------------------------------

outRoot2 = [depthPath, '_layer2'];
outRoot3 = [depthPath, '_layer3'];
outRoot4 = [depthPath, '_layer4'];
outRoot5 = [depthPath, '_layer5'];
outRoot6 = [depthPath, '_layer6'];
outRoot7 = [depthPath, '_layer7'];
outRoot8 = [depthPath, '_layer8'];
    

 
% --------input file names-----------------------------------------------------

dsN = num2str(dataSetNumber);
nCl = num2str(nClusters);
aL = '3'; % abstraction level

vocabulary1Layer = [root, 'statistics/statistics_1_', dsN, '_', nCl, '.mat'];
 
 % files for part selection results
parts3Layer = [root, 'statistics/partsSelectionResults_3_', dsN, '_', nCl, '_a', aL, '.mat'];
parts4Layer = [root, 'statistics/partsSelectionResults_4_', dsN, '_', nCl, '_a', aL, '.mat'];
parts5Layer = [root, 'statistics/partsSelectionResults_5_', dsN, '_', nCl, '_a', aL, '.mat'];
parts6Layer = [root, 'statistics/partsSelectionResults_6_', dsN, '_', nCl, '_a', aL, '.mat'];
parts7Layer = [root, 'statistics/partsSelectionResults_7_', dsN, '_', nCl, '_a', aL, '.mat'];
parts8Layer = [root, 'statistics/partsSelectionResults_8_', dsN, '_', nCl, '_a', aL, '.mat'];

% files for the sieved statistics
statistics1LayerSieved = [root, 'statistics/statisticsSieved_1_', dsN, '_', nCl, '.mat'];
statistics3LayerSieved = [root, 'statistics/statisticsSieved_3_', dsN, '_', nCl, '.mat'];
statistics4LayerSieved = [root, 'statistics/statisticsSieved_4_', dsN, '_', nCl, '.mat'];
statistics5LayerSieved = [root, 'statistics/statisticsSieved_5_', dsN, '_', nCl, '.mat'];
statistics6LayerSieved = [root, 'statistics/statisticsSieved_6_', dsN, '_', nCl, '.mat'];
statistics7LayerSieved = [root, 'statistics/statisticsSieved_7_', dsN, '_', nCl, '.mat'];
statistics8LayerSieved = [root, 'statistics/statisticsSieved_8_', dsN, '_', nCl, '.mat'];

% files for the Aggregated statistics
statistics1LayerAggregated = [root, 'statistics/statisticsAggregated_1_', dsN, '_', nCl, '.mat'];
statistics3LayerAggregated = [root, 'statistics/statisticsAggregated_3_', dsN, '_', nCl, '.mat'];
statistics4LayerAggregated = [root, 'statistics/statisticsAggregated_4_', dsN, '_', nCl, '.mat'];
statistics5LayerAggregated = [root, 'statistics/statisticsAggregated_5_', dsN, '_', nCl, '.mat'];
statistics6LayerAggregated = [root, 'statistics/statisticsAggregated_6_', dsN, '_', nCl, '.mat'];
statistics7LayerAggregated = [root, 'statistics/statisticsAggregated_7_', dsN, '_', nCl, '.mat'];
statistics8LayerAggregated = [root, 'statistics/statisticsAggregated_8_', dsN, '_', nCl, '.mat'];

% files for layer visualization
fileForVisualization3Layer = [root, 'statistics/fileForVisualization_3_', dsN, '_', nCl, '.mat'];
fileForVisualization4Layer = [root, 'statistics/fileForVisualization_4_', dsN, '_', nCl, '.mat'];
fileForVisualization5Layer = [root, 'statistics/fileForVisualization_5_', dsN, '_', nCl, '.mat'];
fileForVisualization6Layer = [root, 'statistics/fileForVisualization_6_', dsN, '_', nCl, '.mat'];
fileForVisualization7Layer = [root, 'statistics/fileForVisualization_7_', dsN, '_', nCl, '.mat'];
fileForVisualization8Layer = [root, 'statistics/fileForVisualization_8_', dsN, '_', nCl, '.mat'];
 
% %--------------------------------------------------------------------------
% % creating a filelist here

if dataSetNumber == 1 || dataSetNumber == 3
    [list_depth, lenF] = extractFileList(fileListPrecomputed, depthPath, depthPathDefault, is_subset, subset_len);
    list_mask = [];
elseif dataSetNumber == 2
    [list_depth, list_mask, ~, lenF] = extractFileListWashington(fileListPrecomputed, depthPath, depthPathDefault, is_subset, subsetPercent);
end


load(vocabulary1Layer); %  cluster1Centres, cluster1Bounds, thresh
depthStep = thresh/4;
%--------------------------------------------------------------------------

if is_first_layer || is_second_layer 
    
    % read the first layer vocabulary

    
    performInference2(list_depth, list_mask, lenF, sigma, sigmaKernelSize, dxKernel, isErrosion, discRadius, nClusters, ...
                areLayersRequired, outRoot2, combs, largestLine, wCoverage, wOverlap, isInhibitionRequired, ...
                is_downsampling, dowsample_rate, dataSetNumber, depthPath, cluster1Centres, cluster1Bounds, thresh, ...
                is_guided, r_guided, eps, is_mask_extended, maxExtThresh1, maxExtThresh2, isReconstructionErrorRequired);
end

% LIST_MASK - TO DELETE
            
if is_third_layer 
    
    elPath = outRoot2; % second layer elements
    is_subset = false;
    
    if dataSetNumber == 2  % extract images with previous layer
        [list_El, ~, ~, lenF] = extractFileListWashingtonForClassification(elPath, is_subset, 1.0);
    else
        [list_El, lenF] = extractFileList(fileListPrecomputed, elPath, '', is_subset, subset_len);
    end
    
    load(statistics3LayerAggregated);  % 'X' ,'frequencies', 'curTS', 'triples'
    load(parts3Layer);   % 'triples3Out', 'coverageOut', 'n3Clusters', 'abstractionLevel');
    load(statistics3LayerSieved);     %   'statistics', 'cluster3Depths', 'outputCoords'
    clear('statistics', 'outputCoords', 'triples', 'frequencies');

    triples3Out = triples3Out(1:n3Clusters, :);
    displacement34 = 6;
    
    performInference3(list_depth, list_El, list_mask, lenF, sigma, sigmaKernelSize, dxKernel, nClusters, n2Clusters,...
                n3Clusters, X, triples3Out, coverageOut, displacement34, ...
                outRoot3, combs, largestLine, wCoverage, wOverlap, isInhibitionRequired, ...
                is_downsampling, dowsample_rate, elPath, meargeThresh3, isErrosion, discRadius, is_guided, r_guided, eps, ...
                is_mask_extended, maxExtThresh1, maxExtThresh2, depthStep, cluster3Depths, dataSetNumber);         
            
end



if is_4th_layer 
    
    elPath = outRoot3; % third layer elements
    
    is_subset = false;
    
    if dataSetNumber == 2  % extract images with previous layer
        [list_El, ~, ~, lenF] = extractFileListWashingtonForClassification(elPath, is_subset, 1.0);
    else
        [list_El, lenF] = extractFileList(fileListPrecomputed, elPath, depthPathDefault, is_subset, subset_len)
    end
    
    load(statistics4LayerAggregated);  % 'X' ,'frequencies', 'curTS', 'triples'
    load(parts3Layer);   % 'triples3Out', 'coverageOut', 'n3Clusters', 'abstractionLevel');
    load(parts4Layer);   % 'triples4Out', 'coverageOut', 'n4Clusters', 'abstractionLevel');
    load(statistics4LayerSieved);     %   'statistics', 'cluster4Depths', 'outputCoords'
    clear('statistics', 'outputCoords', 'triples', 'frequencies', 'triples3Out');

    triples4Out = triples4Out(1:n4Clusters, :);
    
    displacement34 = 6;
    
    performInference4(list_depth, list_El, list_mask, lenF, sigma, sigmaKernelSize, dxKernel, ...
                n3Clusters, n4Clusters, X, triples4Out, coverageOut, displacement34, ...
                outRoot4, combs, largestLine, wCoverage, wOverlap, isInhibitionRequired, ...
                is_downsampling, dowsample_rate, elPath, meargeThresh4, isErrosion, discRadius, is_guided, r_guided, eps, ...
                is_mask_extended, maxExtThresh1, maxExtThresh2, depthStep, cluster4Depths, fileForVisualization3Layer, dataSetNumber);   
       
end



if is_5th_layer 
    
    elPath = outRoot4; % third layer elements
    layerID = 5;
    is_subset = false;
    
    if dataSetNumber == 2  % extract images with previous layer
        [list_El, ~, ~, lenF] = extractFileListWashingtonForClassification(elPath, is_subset, 1.0);
    else
        [list_El, lenF] = extractFileList(fileListPrecomputed, elPath, depthPathDefault, is_subset, subset_len);
    end
    
    load(statistics5LayerAggregated);  % 'X' ,'frequencies', 'curTS', 'triples'
    load(parts4Layer);   % 'triples4Out', 'coverageOut', 'n4Clusters', 'abstractionLevel');
    load(parts5Layer);   % 'triples5Out', 'coverageOut', 'n5Clusters', 'abstractionLevel');
    
    load(statistics5LayerSieved);     %   'statistics', 'cluster5Depths', 'outputCoords'
    clear('statistics', 'outputCoords', 'triples', 'frequencies', 'triples4Out');

    triples5Out = triples5Out(1:n5Clusters, :);
    
    displ3 = 6;
    displ5 = 18;
    displ7 = 54;
    
    performInferenceNext(list_depth, list_El, list_mask, lenF, sigma, sigmaKernelSize, dxKernel, ...
                n4Clusters, n5Clusters, X, triples5Out, coverageOut, displ3, displ5, displ7,  ...
                areLayersRequired, outRoot5, combs, largestLine, wCoverage, wOverlap, isInhibitionRequired, ...
                is_downsampling, dowsample_rate, elPath, meargeThresh5, isErrosion, discRadius, is_guided, r_guided, eps, ...
                is_mask_extended, maxExtThresh1, maxExtThresh2, depthStep, cluster5Depths, fileForVisualization4Layer, ...
                dataSetNumber, layerID);
            
    
end


if is_6th_layer 
    
    elPath = outRoot5; % third layer elements
    layerID = 6;
    is_subset = false;
    
    if dataSetNumber == 2  % extract images with previous layer
        [list_El, ~, ~, lenF] = extractFileListWashingtonForClassification(elPath, is_subset, 1.0);
    else
        [list_El, lenF] = extractFileList(fileListPrecomputed, elPath, depthPathDefault, is_subset, subset_len);
    end
    
    load(statistics6LayerAggregated);  % 'X' ,'frequencies', 'curTS', 'triples'
    load(parts5Layer);   % 'triples5Out', 'coverageOut', 'n5Clusters');
    load(parts6Layer);   % 'triples6Out', 'coverageOut', 'n6Clusters');
    
    load(statistics6LayerSieved);     %   'statistics', 'cluster5Depths', 'outputCoords'
    clear('statistics', 'outputCoords', 'triples', 'frequencies', 'triples5Out');

    triples6Out = triples6Out(1:n6Clusters, :);
    
    displ3 = 6;
    displ5 = 18;
    displ7 = 54;
    
    performInferenceNext(list_depth, list_El, list_mask, lenF, sigma, sigmaKernelSize, dxKernel, ...
                n5Clusters, n6Clusters, X, triples6Out, coverageOut, displ3, displ5, displ7,  ...
                areLayersRequired, outRoot6, combs, largestLine, wCoverage, wOverlap, isInhibitionRequired, ...
                is_downsampling, dowsample_rate, elPath, meargeThresh6, isErrosion, discRadius, is_guided, r_guided, eps, ...
                is_mask_extended, maxExtThresh1, maxExtThresh2, depthStep, cluster6Depths, fileForVisualization5Layer, ...
                dataSetNumber, layerID);
            
    
end




