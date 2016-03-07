% This is a script to learn the 3D compositional hierarchical shape
% vocabulary

% dataSetNumber 
% for Aim@Shape dataSetNumber = 1;list_input
% for Washington dataSetNumber = 2;
% for Vladislav's_STD dataSetNumber = 3;
% This is the version LIBHOP3D_DG - with diffrential geometry


% inputDataType = 1; % depth images
% inputDataType = 2; % meshes
% inputDataType = 3; % point clouds

% partSelectionMethod
% for MDL-based for depth images:              partSelectionMethod = 1 
% for MDL-based for meshes and point clouds:   partSelectionMethod = 2


function [nNClusters] = learnHierarchy()

commonRoot = 'E:/3D/';
root = [commonRoot, 'LibHoP3D_DG/'];
addPaths(root);

dataSetNumber = 5;
dataSetNames = getDataSetNames();

if dataSetNumber < 5
    inputDataType = 1;      % depth images
elseif dataSetNumber == 5
    inputDataType = 2;  % meshes 
else
    inputDataType = 3;  % point clouds
end

nClusters = 1;
n2Clusters = 1;

is_GPU_USED = true;
% if is_GPU_USED
%     g = gpuDevice(1);
%     reset(g);
% end

is_multyScale = true;

if inputDataType == 1
    [iterations, lenSelected, fieldSize, maxRelDepth, quant, partsCoverArea, numSimilar, threshPair, sieve_thresh] = loadPartSelectionParameters(dataSetNumber);
else
    [iterations, lenSelected, fieldSize, maxRelDepth, quant, partsCoverArea, numSimilar, threshPair, sieve_thresh] = loadPartSelectionParametersMesh(dataSetNumber);
end

[statMapPropertiesAll, pairClusteringOptionsAll] = loadStatMapProperties(dataSetNumber);  
[inputPath] = getPathToData(dataSetNumber, commonRoot);
[scales, lineAdders]  = getScales(dataSetNumber, is_multyScale); 
%computeAllscales(scales, inputPath{1,1}, lineAdders, dataSetNumber, inputDataType); % convert all input images to different scales
[receptiveField, offsetsConventional] = getReceptiveFieldSize(dataSetNumber);


%% define output file names

dsN = num2str(dataSetNumber);
nCl = num2str(nClusters);

[partsLayer, fileForVisualizationLayer, partsLayerLoc, partsLayerEnt, partsSpecialSelected, partsLayerAll, calibrationFile] = getStandardFilePaths(root, dsN, nCl);

vocabulary1Layer = [root, 'statistics/statistics_1_', dsN, '_', nCl, '.mat'];
diskFolderVisVoc = [commonRoot, 'Visualized vocabulary/'];
diskFolderStatMap = [commonRoot, 'Visualized statistics/'];

for i = 3:8
    statisticalMapLayer{i}   = [root, 'statistics/statisticalMap_', num2str(i), '_', dsN, '_', nCl, '.mat'];
    folderForLayer{i}        = [diskFolderVisVoc,  dataSetNames{dataSetNumber},'/', num2str(i), '_layer/'];
    folderForLayerStatMap{i} = [diskFolderStatMap, dataSetNames{dataSetNumber},'/', num2str(i), '_layer/'];
end


% if ~exist(vocabulary1Layer, 'file');
%     [cluster1Centres, cluster1Bounds, thresh] = defineFirstLayerQuantiles(nClusters, dataSetNumber);
%     save(vocabulary1Layer, 'cluster1Centres', 'cluster1Bounds', 'thresh');
% else
%     load(vocabulary1Layer); % cluster1Centres, cluster1Bounds, thresh
% end

% [ is_ok ] = VisualizeCornerParts(cluster1Centres); 
% a = 2;

%% define all parameters here
if inputDataType == 1   % depth images
%     [filtOptions, zeroThresh] = loadFilteringParameters(dataSetNumber);
%     options = GetOptions(filtOptions);
end


%% Define what we have to learn
%                                    {1,2,3,4,5,6,7,8,9, 10}
is_layer =                           {0,0,0,0,0,1,1,1,1,1};
is_statisticalMap =                  {0,0,0,0,0,0,0,0,0,0};
is_statistics_collection =           {0,0,0,0,0,1,1,1,1,1};

%                                    {1,2,3,4,5,6,7,8}
is_statistics_sieve_aggregate_Weak = {0,0,0,0,0,1,1,1};
is_localization =                    {0,0,0,0,0,0,0,0};
is_Entropy =                         {0,0,0,0,0,0,0,0};
is_partSelectionSpecialNeeded =      {0,0,0,0,0,0,0,0};

%                                    {1,2,3,4,5,6,7,8}
is_statistics_sieve_aggregate_Strong={0,0,0,0,0,0,0,0};
is_partSelectionNeeded =             {0,0,0,0,0,0,0,0};
combinePartSelection =               {0,0,0,0,1,1,0,0};
visualizeLayer =                     {0,0,0,0,1,1,1,1};
is_inferenceNeeded =                 {0,0,0,0,0,0,0,0};
isDownsampling  =                    {0,0,0,0,0,0,0,0};

%% downsampling of input images

is_downsampling = false;
if is_downsampling && inputDataType == 1
    dowsample_rate = 1.0;
    is_subsetDS = false; % whether we shall use all files for learning
    subsetPercentDS = 1.0;
    UpsampleAllImages(dataSetNumber, inputPath{1,1}, dowsample_rate, is_subsetDS, subsetPercentDS);
end


%% adding zero boundaries

% boundSize = 10;
% parfor i = 1:lenF
%     I = imread(list_input{i});
%     I = I(:,:,1);
%     I = addZeroBoundaries(I, boundSize);
%     imwrite(I, list_input{i}, 'png');
% end



%% inference of the first layer parts in the data

nNClusters{1} = nClusters;
nNClusters{2} = n2Clusters;
nNClusters{3} = 6;
nNClusters{4} = 10 + 1; %% !!!! planar patch
nNClusters{5} = 32; 
% nNClusters{6} = 14 + 1; %% !!!! planar patch
% nNClusters{7} = 13;
% nNClusters{8} = 

for i = 3:8
    tripleOutDepth{i}   = [];
    triplesCurOut{i}    = [];
    partsOutSpecial{i}   = [];
    nNClustersSpecial{i} = [];
end


depthStep = 0.0005;

% maxDist = 2;
% borderOffset34 = 9;
% thresh = 4;

%% calibrate the dataset using an image of the sphere-like objects

if inputDataType == 1  % depth images
    % perform calibration of the dataSet
    
    if exist(calibrationFile, 'file');
        zScale = load(calibrationFile);
        zScale = zScale.zScale;
    else
        % create this file
        zScale = calibrateDataSet(dataSetNumber, inputPath{1,1});
        save(calibrationFile, 'zScale');
    end
end

is_subset = false;
subsetPercent = 1.0;


%% learn the FIRST LEVEL PARTs (non-planar)
% 
if is_layer{1}
    layerID = 1;
    methodID_FL = 1;
    for i = 1:length(scales)
        % get path to the input data at this scale
        str = inputPath{layerID, 1};
        input_path = getPathScale(str, lineAdders{i});
        [list_input, ~, ~, lenFiles] = extractFileListGeneral(input_path, is_subset, subsetPercent, dataSetNumber);
        learnNonPlanarParts(list_input, lenFiles, dataSetNumber, patchRad, methodID_FL);
    end
end

%%

for layerID = 2:8
    
 if is_layer{layerID}
    
    nPrevClusters = nNClusters{layerID-1};
    if is_inferenceNeeded{layerID - 1} 
        infArray = zeros(1, 8);
        infArray(layerID-1) = 1; % inference of the previous layer
        LayersInference(infArray, dataSetNumber, inputDataType, nClusters, is_multyScale);
    end
    
    if layerID == 2
        continue;  % only inference is needed for this layer
    end
              
    str = ['Learning of the layer ', num2str(layerID), '...'];
    disp(str);     
    
    %% collect statistics

    if is_statistics_collection{layerID} 
        
        if inputDataType == 1
            
            % statistics collection for all scales
            for i = 1:length(scales)
                % get path to the input data at this scale
                str = inputPath{layerID, 1};
                input_path = getPathScale(str, lineAdders{i});
                elPath = getElPath(input_path, layerID-1);

                % get a list of elements and depth images
                [list_input, list_mask, ~, lenF] = extractFileListGeneral(input_path, is_subset, subsetPercent, dataSetNumber);
                [list_els] = makeElList(list_input, input_path, elPath);
                [statisticsLayer, ~, ~, ~, ~] = GetStatisticsFiles(root, dsN, nCl, lineAdders{i}, layerID);

                disp('ERROR:TODO');
            end
        end
        
        list_input = {};
        list_els = {};
        
        if inputDataType == 2
                
            % statistics collection for all scales
            for is = 1:length(scales)
                str = inputPath{1, 1};
                input_path = getPathScale(str, lineAdders{is});
                
                [list_input, ~, ~, lenF] = extractFileListGeneral(input_path, is_subset, subsetPercent, dataSetNumber);
                strI21 = inputPath{2,1};
                str2 = getPathScale(strI21, lineAdders{is});
                elPath = getElPath(str2, layerID-1);
                
                list_els = makeElList(list_input, input_path, elPath);
%                 CollectStats_NextLayersMesh_PC(list_els, list_input, lenF, nPrevClusters, layerID, ...
%                                                                     dataSetNumber, receptiveField{layerID}, offsetsConventional{layerID}, statMapPropertiesAll{layerID});
            end

        sieveThresh1{3} = 50; sieveThresh1{4} = 200; sieveThresh1{5} = 20; sieveThresh1{6} = 30; 
        sieveThresh1{7} = 3000;  sieveThresh1{8} = 500; sieveThresh1{9} = 300; sieveThresh1{10}= 300;

    numParts = zeros(1,15);
    for nn = 1:15 
        a{1} = [33, 34, 35, 36, 37];      % box
        a{2} = [1,32,38,40,43,51,57,72];  % bottles
        a{3} = [8 16,30,31];              % bowls
        a{4} = [41, 42];                  % cup
        a{5} = [45, 46, 47, 48, 49, 50, 52];  % cutting board
        a{6} = [39];                      % can
        a{7} = [13, 14 15 44];            % tray
        a{8} = [2 3 4 73 74 75 76];       % plate 
        a{9} = [6,7,9];  % tea cup
        a{10} = [5];  % salt
        a{11} = [10, 11 12, 65];  % teapot
        a{12} = [17:29 60];  % vase
        a{13} = [53:56 58 59]; % frying pan
        a{14} = [61:64, 66]; % Jug
        a{15} = [67:71]; % mug

        aa = [a{1:nn}];
        lenF = length(aa);

        if nn == 1
            prev = 1;
        else 
            prev = length([a{1:nn-1}]) + 1;
        end


        [statMapLeft, statMapRight] = buildStatMap_NextLayers(aa, prev, lenF, nPrevClusters, offsetsConventional{layerID}, statMapPropertiesAll{layerID}, layerID,...
                                                                lenF * sieveThresh1{layerID}, is_GPU_USED);                                                   
         save('Temp/statMap.mat', 'statMapLeft', 'statMapRight', '-v7.3');

%               clustering of the staistical maps
        [pairsLeft, pairsRight] = statMapClustering(nPrevClusters, statMapPropertiesAll{layerID}, pairClusteringOptionsAll{layerID}, offsetsConventional{layerID}, layerID, lenF);
        save('Temp/pairsLR.mat', 'pairsLeft', 'pairsRight');

        numParts(nn) = size(pairsLeft, 1) + size(pairsRight, 1);

    end

    dd = load('Temp/pairsLR.mat'); 
    pairsLeft = dd.pairsLeft;
    pairsRight = dd.pairsRight;
            
%              % assign each point in the file with statistics by the ID

%             tic
%      assignIdsToStatisticsFile(pairsLeft, pairsRight, pairClusteringOptionsAll{layerID}, nPrevClusters, lenF, statMapPropertiesAll{layerID}, layerID, is_GPU_USED);
% %             toc
%             
%              % aggregate statistics of the layers
%             nPrevClusters = size(pairsLeft, 1) + size(pairsRight, 1);
%             [X, frequencies, triples] = Aggregate_statVI(lenF, nPrevClusters, sieveThresh1{layerID}, pairsLeft, pairsRight);
%              
%             % perform part selection at this layer
%             iterations = 40;
%             lenSelected{3} = 100; lenSelected{4} = 100; lenSelected{5} = 250; lenSelected{6} = 230; lenSelected{7} = 250; lenSelected{8} = 250;
%             lenSelected{9} = 150; lenSelected{10} = 150;
%             
%              
%             ids = frequencies >= log(lenF) * lenSelected{layerID};
%             frequencies = frequencies(ids);
%             X = X(ids, :);
% %             D = computePartPairWiseDistances(X, layerID, pairsLeft, pairsRight, receptiveField{2}, true);
% 
%             aa = load('Temp/pairwiseDist.mat');
%             D = aa.D;
% 
%             [partsOut, lenOut, ORTable] = PartSelectionMeshFreq_VI(X, frequencies, D, layerID);
% % % %              [partsOut, coverageOut, lenOut, ORTable] = PartSelectionMeshMDL_VI(list_els, lenF, nPrevClusters, iterations, lenSelected{layerID}, D, layerID);
%             
%             partsOut = partsOut(1:lenOut, :);
% % % % % % %             coverageOut = coverageOut(1:lenOut);
%             pairsAll = [pairsLeft; pairsRight];
%             nNClusters{layerID} = lenOut;
%             if layerID == 4 || layerID == 6  % include a planar part of the next scale
%                 nNClusters{layerID} = nNClusters{layerID} + 1;
%             end    
%             strOut = ['Temp/','layer', num2str(layerID), '/partSelection', num2str(layerID), '.mat'];
%             save(strOut, 'partsOut', 'pairsAll', 'nNClusters');
%             strLayer = ['Temp/OR_node_layer_', num2str(layerID), '.mat' ];
%             save(strLayer, 'ORTable');
%             
%             ReadVocabulary(layerID);
%             if layerID == 4 || layerID == 6
%                 VisualizeLayerVI(layerID, 2); % starting from element 2
%             else
                VisualizeLayerVI(layerID);
%             end
% % %             inference of the next layer based on statistics
            MakeInfrerenceVICorrect(list_els, list_input, lenF, nPrevClusters, layerID, ...
                                                dataSetNumber, receptiveField{layerID}, offsetsConventional{layerID}, statMapPropertiesAll{layerID});

% % % % %                layerInferenceNext_VI(list_input, list_els, lenF, layerID, is_GPU_USED);
             a = 2;
        end
                
%             save(statisticsLayer, 'outputStatistics', 'outputCoords', 'outputFrames', '-v7.3');
%             save('Temp/depth_files.mat', 'list_input', 'list_mask', 'list_els');
%             a = 2;
    end
    
    load('Temp/depth_files.mat');
    
    %% Weak sieve and aggregation of the statistics

    if is_statistics_sieve_aggregate_Weak{layerID}
        
        weak_multiplier = 0.01;
        
        for i = 1:length(scales)
            
            % load statistics
            [statisticsLayer, statisticsLayerSieved_Weak, statisticsLayerAggregated_Weak, ~, ~] = GetStatisticsFiles(root, dsN, nCl, lineAdders{i}, layerID);
            load(statisticsLayer);
            
            [outputStatistics, ~] = DiscretizeOrientations(outputStatistics, layerID, curTS, nClusters);

            % sieve
            numDisps = 2;
            pairThresh = weak_multiplier *curTS * threshPair{layerID};
            [statistics, outputCoords, outputFrames, clusterCurDepths, curTS] = sieveStatistics(outputStatistics, outputCoords, outputFrames, numDisps, ...
                                                            pairThresh, nPrevClusters, quant{layerID}, maxRelDepth{layerID}, is_GPU_USED);
            clear('outputStatistics');

            % aggregate
            sieveThresh = weak_multiplier *curTS * sieve_thresh{layerID};
            is_GPU_USED = false;
            [X, frequencies, triples, is_sparse] = Aggregate_3Layer(statistics, nPrevClusters, curTS, sieveThresh, is_GPU_USED);
            [ind, statistics] = Sieve3afterAggregation(statistics, triples, sieveThresh, is_sparse, is_GPU_USED);
            outputCoords = outputCoords(ind, :);
            outputFrames = outputFrames(ind, :);
            curTS = size(outputCoords, 1);

            % save the processed statistics

            save(statisticsLayerSieved_Weak, 'statistics', 'clusterCurDepths', 'outputCoords', 'outputFrames', '-v7.3');
            save(statisticsLayerAggregated_Weak, 'X' ,'frequencies', 'curTS', 'triples', '-v7.3');
        end
    end
    
    
  %% compute additional part selection criteria  
  
    if is_Entropy{layerID} % part selection based on discrimination criteria
        
        [partEntropy] = computePartsEntropy(list_input, statisticsLayerSieved_Weak{layerID}, statisticsLayerAggregated_Weak{layerID}, nPrevClusters+1, dataSetNumber);
        save('Temp/entropy.mat', 'partEntropy');
        
        load(statisticsLayerAggregated_Weak{layerID}); 
        
%     Next lines come for visualization
%         [partEntropy, inds] = sort(partEntropy, 'ascend');
%         X = X(inds,:);
%         a = [partEntropy, X];
        
        inds = partEntropy < 1.50;
        triplesOutEnt = X(inds,:);
        partEntropyOut = partEntropy(inds);
        numSelectedEnt = size(triplesOutEnt, 1);
        save(partsLayerEnt{layerID}, 'triplesOutEnt', 'numSelectedEnt', 'partEntropyOut');
    end
        
    if is_localization{layerID} % part selection based on localization criteria
        
%         if ~exist('partEntropy', 'var')
%             load('Temp/entropy.mat');
%         end
        
        winSize = 6;
        [localizationOut] = computePartsLocalization(list_input, statisticsLayerSieved_Weak{layerID}, statisticsLayerAggregated_Weak{layerID}, nPrevClusters+1, winSize);
%         save('Temp/localization.mat', 'localizationOut', 'partEntropy');           
        load(statisticsLayerAggregated_Weak{layerID});
        
%     Next lines come for visualization
        [localizationOut, inds] = sort(localizationOut, 'ascend');
        X = X(inds,:);
        a = [localizationOut,X];

        inds = localizationOut >= 0.25;
        triplesOutLoc = X(inds,:);
        numSelectedLoc = size(triplesOutLoc, 1);
        
        partLocEntropyOut = 5*ones(numSelectedLoc, 1);
        save(partsLayerLoc{layerID}, 'triplesOutLoc', 'numSelectedLoc', 'partLocEntropyOut');
        
    end
    
    
    %% localization-discriminative-based part selection procedure  
    
    if is_partSelectionSpecialNeeded{layerID}
        if ~exist('clusterCurDepths')
            bb = load(statisticsLayerSieved_Weak{layerID});
            clusterCurDepths = bb.clusterCurDepths;
        end
        [partsOutSpecial{layerID}, nNClustersSpecial{layerID}, partEntropySpecial] = partSelectionSpecial(false, true, partsLayerEnt{layerID},...
            partsLayerLoc{layerID},  nClusters, fileForVisualizationLayer, clusterCurDepths, displ, ...
            layerID, cluster1Centres, depthStep, partsCoverArea{layerID});          
  
        save(partsSpecialSelected{layerID}, 'partsOutSpecial', 'nNClustersSpecial', 'partEntropySpecial');
    end
    

    %% Strong Sieve and Aggregation of the statistics
    
    if is_statistics_sieve_aggregate_Strong{layerID} 
        
        strongMultiplier = 0.05;
        
        for i = 1:length(scales)
            
            [statisticsLayer, ~, ~, statisticsLayerSieved, statisticsLayerAggregated] = GetStatisticsFiles(root, dsN, nCl, lineAdders{i}, layerID);
            load(statisticsLayer);

            [outputStatistics, ~] = DiscretizeOrientations(outputStatistics, layerID, curTS, nClusters);
            % sieve
            pairThresh = strongMultiplier * curTS * threshPair{layerID};
            numDisps = 2;

            [statistics, outputCoords, outputFrames, clusterCurDepths, curTS] = sieveStatistics(outputStatistics, outputCoords, outputFrames, numDisps, ...
                                                pairThresh, nPrevClusters, quant{layerID}, maxRelDepth{layerID}, is_GPU_USED);
            clear('outputStatistics');

            % aggregate        
            sieveThresh = strongMultiplier * curTS * sieve_thresh{layerID};
            [X, frequencies, triples, is_sparse] = Aggregate_3Layer(statistics, nPrevClusters, curTS, sieveThresh, is_GPU_USED);
            [ind, statistics] = Sieve3afterAggregation(statistics, triples, sieveThresh, is_sparse, is_GPU_USED);
            outputCoords = outputCoords(ind, :);
            outputFrames = outputFrames(ind, :);
            curTS = size(outputCoords, 1);

            save(statisticsLayerSieved, 'statistics', 'clusterCurDepths', 'outputCoords', 'outputFrames', '-v7.3');
            save(statisticsLayerAggregated, 'X' ,'frequencies', 'curTS', 'triples', '-v7.3');
        end
        
    end 
    
    
    
 %% Coverage-based part selection procedure   
         
    if is_partSelectionNeeded{layerID}
            
        if inputDataType == 1  % depth images
            
            [triplesCurOut{layerID}, coverageOut, nNClusters{layerID}] = PartSelectionFull(nClusters, nPrevClusters, dataSetNumber, partsCoverArea{layerID}, ... 
                iterations{layerID}, layerID, fileForVisualizationLayer, lenSelected{layerID}, cluster1Centres, numSimilar{layerID}, offsetsConventional, depthStep, ...
                is_multyScale, scales, lineAdders, inputPath{layerID, 1}, is_subset, subsetPercent, root, dsN, nCl);
            
        elseif inputDataType == 2  % meshes
            for i = 1:length(scales)
                % get path to the input data at this scale
                str = inputPath{layerID, 1};
                input_path = getPathScale(str, lineAdders{i});
                elPath = getElPath(input_path, layerID-1);

                % get a list of elements and depth images
                [list_input, list_mask, ~, lenF] = extractFileListGeneral(input_path, is_subset, subsetPercent, dataSetNumber);
                [list_els] = makeElList(list_input, input_path, elPath);

                [statisticsLayer, ~, ~, statisticsLayerSieved, statisticsLayerAggregated] = GetStatisticsFiles(root, dsN, nCl, lineAdders{i}, layerID);
                [triplesCurOut{layerID}, coverageOut,nNClusters{layerID}] = PartSelectionMeshMDL(statisticsLayer, statisticsLayerSieved, statisticsLayerAggregated, list_input, list_els, layerID,...
                            nClusters, n2Clusters, fileForVisualizationLayer, offsetsConventional, depthStep, fieldSize{layerID}, cluster1Centres, iterations{layerID}, ...
                            numSimilar{layerID}, lenSelected{layerID});
            end
        end

        partEntropyMain = 5*ones(nNClusters{layerID},1);    
        % previously selected 
        save(partsLayer{layerID}, 'triplesCurOut', 'coverageOut', 'nNClusters', 'partEntropyMain');            
                  

    end
    
    % store the visualization of the vocabulary to the folder
    
    if combinePartSelection{layerID}
        
        try
            load(partsSpecialSelected{layerID});       % 'partsOutSpecial', 'nNClustersSpecial', 'partEntropySpecial'
        
            % this is to compute part's entropy for triplesCurOut
            partsEntropy{layerID} = [partEntropySpecial];
            triplesCurOut{layerID} = [partsOutSpecial{layerID}];
            nNClusters{layerID} = nNClustersSpecial{layerID};

            save(partsLayerAll{layerID}, 'triplesCurOut', 'nNClusters', 'partsEntropy');
        catch exception
    
            % if there on of the files does not exist
            load(partsLayer{layerID});
            partsEntropy = partEntropyMain;
            save(partsLayerAll{layerID}, 'triplesCurOut', 'nNClusters', 'partsEntropy');
        end
    end
    
    
    %% visualization of parts of the layer
    
    if visualizeLayer{layerID}

        load(partsLayerAll{layerID});
        [~, ~, ~, statisticsLayerSieved, ~] = GetStatisticsFiles(root, dsN, nCl, lineAdders{1}, layerID);
        load(statisticsLayerSieved);      %  'statistics', 'clusterCurDepths', 'outputCoords'
        
        if mod(layerID,2) == 1
            LI = layerID + 1;  % to make a squared visualization field
        else
            LI = layerID;
        end
            
        % load(fileForVisualization3Layer);

        % prepare data for visualization (convert to the right format)
        if layerID == 3
            tripleOutDepth{layerID} = store3Layer(triplesCurOut{layerID}, clusterCurDepths, nNClusters{layerID}, nClusters);
        else
            tripleOutDepth{layerID} = store4Layer(triplesCurOut{layerID}, clusterCurDepths, nNClusters{layerID}, nClusters);
        end
        
        save(fileForVisualizationLayer{layerID}, 'tripleOutDepth');

        [is_ok] = layer_N_demonstrator(layerID, tripleOutDepth, offsetsConventional, ...
                                            nClusters, nNClusters{layerID}, folderForLayer{layerID}, fieldSize{LI}, depthStep, cluster1Centres, isFIG, receptiveField{2}/2);
%                                         (layerID, tripleOutDepth, offsetsConventional, ...
%                                             nClusters, nNClusters, str_folder, fieldSize, depthStep, cluster1Centres, isFIG, elementRadius)
    end 
  end
end

end




% % %             if layerID >= 3 %|| layerID == 6
% % %                 ORTable = MakeOrNode(layerID, receptiveField{2});
% % %                 strLayer = ['Temp/OR_node_layer_', num2str(layerID), '.mat' ];
% % %                 save(strLayer, 'ORTable');
% % %             end
