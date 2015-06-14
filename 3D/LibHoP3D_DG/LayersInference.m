% this is to perform inference of all levels of the hierarchy

% do not forget: 'matlabpool open 4' for parallel computing

% dataSetNumber 
% for Aim@Shape datasetnumber = 1;
% for Washington datasetnumber = 2

function [] = LayersInference(infArray, dataSetNumber, inputDataType, nClusters, is_multyScale)


if nargin == 0
    % --------Define the layers for inference--------------------------------
    infArray = [1,0,0,0,0,0,0,0];
    dataSetNumber = 5;
    nClusters = 7;
    inputDataType = 2;
end

displ3 = 6;
displ5 = 18;
displ7 = 52;
n2Clusters = nClusters^2;
is_GPU_USED = false;

if is_GPU_USED
    g = gpuDevice(1);
    reset(g);
end

% define folders configuration
commonRoot = 'D:/';
root = [commonRoot, 'LibHoP3D_DG/'];
addPaths(root);

[inputPath] = getPathToData(dataSetNumber, commonRoot);
downSamplingFactor = 3;  % reduction of resolution


is_subset = false; % whether we shall use all files for learning
subsetPercent = 1.0;
is_downsampling = false;
dowsample_rate = 1;


% define all filtering parameters 
receptiveField = getReceptiveFieldSize(dataSetNumber);
if is_multyScale
    [scales, lineAdders] = getScales(dataSetNumber, is_multyScale);
    lenSc = length(scales);
end

if inputDataType == 1 % depth images
    [filtOptions, zeroThresh] = loadFilteringParameters(dataSetNumber);
end


% [~, ~, fieldSize, ~, ~, ~, numSimilar, ~, ~] = loadPartSelectionParameters(dataSetNumber);
% [~, ~, inferenceElType, inferenceElRadius] = loadLearningInferenceStructElement(dataSetNumber);

[ meargeThresh ] = defineMeargeThreshes(1, dataSetNumber);

%                           1,2,3,4,5,6,7,8
is_inhibition    =         {0,0,0,0,0,0,0,0};
is_reconstructionError =   {0,0,0,0,0,0,0,0};
isDownsampling =           {0,0,0,0,0,0,0,0};

% --------input file names-------------------------------------------------

dsN = num2str(dataSetNumber);
nCl = num2str(nClusters);

is_overwrite = true;


vocabulary1Layer = [root, 'statistics/statistics_1_', dsN, '_', nCl, '.mat'];
load(vocabulary1Layer); %  cluster1Centres, cluster1Bounds, thresh
% depthStep = thresh/4;

% [~, statisticsLayerSieved_Strong, statisticsLayerAggregated_Strong, statisticsLayerSieved_Weak, statisticsLayerAggregated_Weak, ...
%     ~, fileForVisualizationLayer, ~, ~, ~, partsLayerAll, calibrationFile] = getStandardFilePaths(root, dsN, nCl);

[~, fileForVisualizationLayer, ~, ~, ~, partsLayerAll, calibrationFile] = getStandardFilePaths(root, dsN, nCl);

zScale = load(calibrationFile);
zScale = zScale.zScale;  % ratio of z and x dimensions


if infArray(2) % perform inference of the second layer 
    
    disp('Inference of the 2nd layer ...');
    
    if is_multyScale  % perform inference at each scale separately
        
        for i = lenSc:-1:1  % For each scale perform inference separately
            str = inputPath{1, 1};
            str1 = getPathScale(str, lineAdders{i});
            [list_input, ~, ~, lenF] = extractFileListGeneral(str1, is_subset, subsetPercent, dataSetNumber);
            
            strI21 = inputPath{2,1};
            str2 = getPathScale(strI21, lineAdders{i});
            strE = getElPath(str2, 2);
            performInference2DG(list_input, lenF, inputDataType, receptiveField{1}, str2, zScale, filtOptions, is_overwrite, strE);
                                           
        end
    end
end


% LIST_MASK - TO DELETE

for layerID = 3:6
    
    if infArray(layerID)
        
        str = ['Inference of the layer ', num2str(layerID), '...'];
        disp(str); 

        [list_input, list_mask, ~, lenF] = extractFileListGeneral(inputPath{layerID-1, 2}, is_subset, subsetPercent, dataSetNumber); 
        [list_els] = makeElList(list_input, inputPath{layerID-1, 2}, elPath{layerID-1, 2});
        
        try   
            load(statisticsLayerAggregated_Weak{layerID});  % 'X' ,'frequencies', 'curTS', 'triples'
            load(statisticsLayerSieved_Weak{layerID});     %   'statistics', 'clusterCurDepths', 'outputCoords'
        catch exception  % if this statistics is not computed
            load(statisticsLayerSieved_Strong{layerID});     %   'statistics', 'clusterCurDepths', 'outputCoords'
            load(statisticsLayerAggregated_Strong{layerID});  % 'X' ,'frequencies', 'curTS', 'triples'
        end
        
        load(partsLayerAll{layerID});   % 'triplesCurOut', 'nNClusters', 'partEntropy'
        clear('statistics', 'outputCoords', 'triples', 'frequencies');

        triplesCurOut{layerID} = triplesCurOut{layerID}(1:nNClusters{layerID}, :);

        performInferenceNext(list_input, list_els, list_mask, lenF, sigma, sigmaKernelSize, dxKernel, ...
                    nNClusters{layerID-1}, nNClusters{layerID}, nClusters, X, triplesCurOut{layerID}, partsEntropy{layerID}, displ3, displ5, displ7, ...
                    elPath{layerID, 1}, is_inhibition, inferenceElType{layerID}, inferenceElRadius{layerID}, ...  
                    is_downsampling, dowsample_rate, elPath{layerID-1, 2}, meargeThresh{layerID}, isErrosion, discRadius, is_guided, r_guided, eps, ...
                    is_mask_extended, maxExtThresh1, maxExtThresh2, depthStep, clusterCurDepths, fileForVisualizationLayer{layerID-1}, ...
                    dataSetNumber, layerID, cluster1Centres, fieldSize{layerID}, cluster1Bounds, numSimilar{layerID}, is_overwrite);

        % count number of detections for each part

        [list_els, ~, ~, lenF] = extractFileListForClassificationGeneral(elPath{layerID, 1}, is_subset, 1.0, dataSetNumber);
        tableEls = ComputeStatsAfterInference(list_els, nNClusters{layerID});
        tableEls = tableEls';
        a = [triplesCurOut{layerID}, tableEls];
        
            % perform downsampling immediately after inference procedure
    
        if isDownsampling{layerID}                                                     % downsamples MARKS images in place!!!

            marksDownsampling(list_els, lenF, downSamplingFactor, elPath{layerID, 1}, elPath{layerID, 2});
            
            % check statistics after downsampling
            [list_els, ~, ~, lenF] = extractFileListForClassificationGeneral(elPath{layerID, 2}, is_subset, 1.0, dataSetNumber);
            [ tableEls ] = ComputeStatsAfterInference(list_els, nNClusters{layerID});
            tableEls = tableEls';
            a = [triplesCurOut{layerID}, tableEls];

        end

    end
end
end

function sP = getPathScale(pathBase, adder)
    sP = [pathBase, '_', adder];
end

function eP = getElPath(in_p, l_ID)
    eP = [in_p, '_layer', num2str(l_ID)];
end


