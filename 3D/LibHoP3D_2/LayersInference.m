% this is to perform inference of all levels of the hierarchy

% do not forget: 'matlabpool open 4' for parallel computing

% dataSetNumber 
% for Aim@Shape datasetnumber = 1;
% for Washington datasetnumber = 2

function [] = LayersInference(infArray, dataSetNumber, nClusters)

if nargin == 0
    % --------Define the layers for inference--------------------------------
    infArray = [0,0,0,0,0,0,0,0];
    dataSetNumber = 3;
    nClusters = 7;
end

displ3 = 6;
displ5 = 18;
displ7 = 52;
isFIG = false;

n2Clusters = nClusters^2;

% % here we initialize the perallel computing
% g = gpuDevice(1);
% reset(g);

% matlabpool open 10    % for parallel computin

% define folders configuration
commonRoot = 'D:/'; 
root = [commonRoot, 'LibHoP3D/'];
addPaths(root);

[depthPath, outRoot ] = getPathToData(dataSetNumber, commonRoot);
downSamplingFactor = 3;  % reduction of resolution

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

[dxKernel, dyKernelTop, dyKernelBottom, dxKernelBack, dxKernelForward, combs, largestLine, sigma, sigmaKernelSize, isErrosion, discRadius, is_guided, r_guided, eps, ...
    is_mask_extended, maxExtThresh1, maxExtThresh2] = loadFilteringParameters(dataSetNumber);

[~, ~, fieldSize, ~, ~, ~, numSimilar, ~, ~] = loadPartSelectionParameters(dataSetNumber);
[~, ~, inferenceElType, inferenceElRadius] = loadLearningInferenceStructElement(dataSetNumber);



% ------------ parameters for line discretization---------------------
% min [error + wOverlap * overlap - wCoverage*coverage ]
if dataSetNumber == 1 || dataSetNumber == 3
    wCoverage = 0.25;
    wOverlap = 0.4;
elseif dataSetNumber == 2    
    wCoverage = 0.05;
    wOverlap = 0.05;
end

[ meargeThresh ] = defineMeargeThreshes(1, dataSetNumber);

isMaskResize = false;

is_inhibition{1} = 0;
is_inhibition{2} = 0; 
is_inhibition{3} = 0;
is_inhibition{4} = 0; 
is_inhibition{5} = 0; 
is_inhibition{6} = 0;
is_inhibition{7} = 0; 
is_inhibition{8} = 0;

is_reconstructionError{1} = false;
is_reconstructionError{2} = false; 
is_reconstructionError{3} = false;
is_reconstructionError{4} = false; 
is_reconstructionError{5} = false; 
is_reconstructionError{6} = false;
is_reconstructionError{7} = false; 
is_reconstructionError{8} = false;

isDownsampling{4} = 0;
isDownsampling{6} = 0;

 
% --------input file names-------------------------------------------------

dsN = num2str(dataSetNumber);
nCl = num2str(nClusters);
aL = '3'; % abstraction level

vocabulary1Layer = [root, 'statistics/statistics_1_', dsN, '_', nCl, '.mat'];

[~, ~, ~, statisticsLayerSieved_Weak, statisticsLayerAggregated_Weak, ...
    ~, fileForVisualizationLayer, ~, ~, ~, partsLayerAll] = getStandardFilePaths(root, dsN, nCl, aL);

load(vocabulary1Layer); %  cluster1Centres, cluster1Bounds, thresh
depthStep = thresh/4;

%--------------------------------------------------------------------------

if infArray(1) || infArray(2) 
    
    % read the first layer vocabulary
    disp('Inference of the 2nd layer ...');  

    if dataSetNumber == 1 || dataSetNumber == 3
        [list_depth, lenF] = extractFileList(depthPath{1}, is_subset, subset_len);
        list_mask = [];
    elseif dataSetNumber == 2
        [list_depth, list_mask, ~, lenF] = extractFileListWashington(depthPath{1}, is_subset, subsetPercent);
    end
    
    performInference2(list_depth, list_mask, lenF, sigma, sigmaKernelSize, dxKernel, isErrosion, discRadius, nClusters, ...
                infArray, outRoot{2}, combs, largestLine, wCoverage, wOverlap, is_inhibition, ...
                is_downsampling, dowsample_rate, dataSetNumber, depthPath{2}, cluster1Centres, cluster1Bounds, thresh, ...
                is_guided, r_guided, eps, is_mask_extended, maxExtThresh1, maxExtThresh2, is_reconstructionError, ...
                dyKernelTop, dyKernelBottom, dxKernelBack, dxKernelForward);
    
            
    % count number of detections for each part
    
    if dataSetNumber == 2  % extract images with previous layer
        [list_els, ~, ~, lenF] = extractFileListWashingtonForClassification(outRoot{2}, is_subset, 1.0);
    else
        [list_els, lenF] = extractFileList(outRoot{2}, is_subset, subset_len);
    end
    
    [ tableEls ] = ComputeStatsAfterInference(list_els, n2Clusters);
    
end

% LIST_MASK - TO DELETE

for layerID = 3:6
    
    if infArray(layerID)
        
        str = ['Inference of the layer ', num2str(layerID), '...'];
        disp(str); 
        
        if dataSetNumber == 1 || dataSetNumber == 3
            [list_depth, lenF] = extractFileList(depthPath{layerID}, is_subset, subset_len);
            list_mask = [];
        elseif dataSetNumber == 2
            [list_depth, list_mask, ~, lenF] = extractFileListWashington(depthPath{layerID}, is_subset, subsetPercent);
        end
         
        [list_els] = makeElList(list_depth, depthPath{layerID-1}, outRoot{layerID-1});

        load(statisticsLayerAggregated_Weak{layerID});  % 'X' ,'frequencies', 'curTS', 'triples'
        load(partsLayerAll{layerID});   % 'triplesCurOut', 'nNClusters', 'partEntropy'
        load(statisticsLayerSieved_Weak{layerID});     %   'statistics', 'clusterCurDepths', 'outputCoords'
        clear('statistics', 'outputCoords', 'triples', 'frequencies');

        triplesCurOut{layerID} = triplesCurOut{layerID}(1:nNClusters{layerID}, :);

        performInferenceNext(list_depth, list_els, list_mask, lenF, sigma, sigmaKernelSize, dxKernel, ...
                    nNClusters{layerID-1}, nNClusters{layerID}, nClusters, X, triplesCurOut{layerID}, partsEntropy, displ3, displ5, displ7, ...
                    outRoot{layerID}, is_inhibition, inferenceElType{layerID}, inferenceElRadius{layerID}, ...  
                    is_downsampling, dowsample_rate, outRoot{layerID-1}, meargeThresh{layerID}, isErrosion, discRadius, is_guided, r_guided, eps, ...
                    is_mask_extended, maxExtThresh1, maxExtThresh2, depthStep, clusterCurDepths, fileForVisualizationLayer{layerID-1}, ...
                    dataSetNumber, layerID, cluster1Centres, fieldSize{layerID}, cluster1Bounds, numSimilar{layerID});

        % count number of detections for each part

        if dataSetNumber == 2  % extract images with previous layer
            [list_els, ~, ~, lenF] = extractFileListWashingtonForClassification(outRoot{layerID}, is_subset, 1.0);
        else
            [list_els, lenF] = extractFileList(outRoot{layerID}, is_subset, subset_len);
        end

        [ tableEls ] = ComputeStatsAfterInference(list_els, nNClusters{layerID});
        tableEls = tableEls';
        a = [triplesCurOut{layerID}, tableEls];
        
            % perform downsampling immediately after inference procedure
    
        if isDownsampling{layerID}                                                     % downsamples MARKS images in place!!!

            marksDownsampling(list_depth, list_els, list_mask, lenF, dataSetNumber, downSamplingFactor, isErrosion, discRadius, ...
                is_mask_extended, maxExtThresh1, maxExtThresh2, isMaskResize);

            [ tableEls ] = ComputeStatsAfterInference(list_els, nNClusters{layerID});
            tableEls = tableEls';
            a = [triplesCurOut{layerID}, tableEls];

        end

    end
end

   
% 
% %% change path to the data (now use a path to downsampled data)
% 
% 
% 
% 
% 
% %% inference of the next layers
% 
% if is_5th_layer 
%     
%     disp('Inference of the 5th layer ...');
%     
%     elPath = outRoot4; % third layer elements
%     layerID = 5;
%     is_subset = false;
%     
%     fieldSize = [53, 17, 251];
%     
%     if dataSetNumber == 2  % extract images with previous layer
%         [list_els, ~, ~, lenF] = extractFileListWashingtonForClassification(elPath, is_subset, 1.0);
%     else
%         [list_els, lenF] = extractFileList(fileListPrecomputed, elPath, depthPathDefault, is_subset, subset_len);
%     end
%     
%     load(statistics5LayerAggregated);  % 'X' ,'frequencies', 'curTS', 'triples'
%     load(parts4Layer);   % 'triples4Out', 'coverageOut', 'n4Clusters', 'abstractionLevel');
%     load(parts5Layer);   % 'triples5Out', 'coverageOut', 'n5Clusters', 'abstractionLevel');
%     
%     load(statistics5LayerSieved);     %   'statistics', 'cluster5Depths', 'outputCoords'
%     clear('statistics', 'outputCoords', 'triples', 'frequencies', 'triples4Out');
% 
%     triples5Out = triples5Out(1:n5Clusters, :);
%     
%     displ3 = 6;
%     displ5 = 18;
%     displ7 = 54;
%     
%     performInferenceNext(list_depth, list_els, list_mask, lenF, sigma, sigmaKernelSize, dxKernel, ...
%                 n4Clusters, n5Clusters, nClusters, X, triples5Out, coverageOut, displ3, displ5, displ7,  ...
%                 outRoot5, combs, largestLine, wCoverage, wOverlap, isInhibitionRequired, ...
%                 is_downsampling, dowsample_rate, elPath, meargeThresh{layerID}, isErrosion, discRadius, is_guided, r_guided, eps, ...
%                 is_mask_extended, maxExtThresh1, maxExtThresh2, depthStep, cluster5Depths, fileForVisualization4Layer, ...
%                 dataSetNumber, layerID, cluster1Centres, fieldSize, thresh);
%             
%             
%     if dataSetNumber == 2  % extract images with layer 4 parts
%         [list_els, ~, ~, lenF] = extractFileListWashingtonForClassification(outRoot5, is_subset, 1.0);
%     else
%         [list_els, lenF] = extractFileList(fileListPrecomputed, outRoot5, depthPathDefault, false, subset_len);
%     end
% 
%     [ tableEls ] = ComputeStatsAfterInference(list_els, n5Clusters);
%     tableEls = tableEls';
%     a = [triples5Out, tableEls];
%      
%     
% end
% 
% 
% if is_6th_layer 
%     
%     disp('Inference of the 6th layer ...');
%     
%     elPath = outRoot5; % third layer elements
%     layerID = 6;
%     is_subset = false;
%     fieldSize = [53, 53, 351];
%     
%     if dataSetNumber == 2  % extract images with previous layer
%         [list_els, ~, ~, lenF] = extractFileListWashingtonForClassification(elPath, is_subset, 1.0);
%     else
%         [list_els, lenF] = extractFileList(fileListPrecomputed, elPath, depthPathDefault, is_subset, subset_len);
%     end
%     
%     load(statistics6LayerAggregated);  % 'X' ,'frequencies', 'curTS', 'triples'
%     load(parts5Layer);   % 'triples5Out', 'coverageOut', 'n5Clusters');
%     load(parts6Layer);   % 'triples6Out', 'coverageOut', 'n6Clusters');
%     
%     load(statistics6LayerSieved);     %   'statistics', 'cluster5Depths', 'outputCoords'
%     clear('statistics', 'outputCoords', 'triples', 'frequencies', 'triples5Out');
% 
%     triples6Out = triples6Out(1:n6Clusters, :);
%     
%     displ3 = 6;
%     displ5 = 18;
%     displ7 = 54;
%     
%     
%     performInferenceNext(list_depth, list_els, list_mask, lenF, sigma, sigmaKernelSize, dxKernel, ...
%                 n5Clusters, n6Clusters, nClusters, X, triples6Out, coverageOut, displ3, displ5, displ7,  ...
%                 outRoot6, combs, largestLine, wCoverage, wOverlap, isInhibitionRequired, ...
%                 is_downsampling, dowsample_rate, elPath, meargeThresh{layerID}, isErrosion, discRadius, is_guided, r_guided, eps, ...
%                 is_mask_extended, maxExtThresh1, maxExtThresh2, depthStep, cluster6Depths, fileForVisualization5Layer, ...
%                 dataSetNumber, layerID, cluster1Centres, fieldSize, thresh);
%             
%     
% end
% 
% 


