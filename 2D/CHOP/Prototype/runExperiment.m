%> Name: runExperiment
%>
%> Description: The entry function to libHOP code. Calls libHOP with
%> desired parameters on specified dataset. For each image in the dataset,
%> we create a graph out of level 1 gabor filter responses, and form
%> edges. Then, we compress samples of each category with SUBDUE to get 
%>
%> @param datasetName Name of the dataset to work on. 
%> @param imageExtension The extension of the files to work on. Examples
%> include '.jpg', '.png', '_crop.png'...
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 18.11.2013
%> Ver 1.1 on 05.12.2013 Various parameter additions, 'mode' changes
%> Ver 1.2 on 12.01.2014 Comment changes for unified code look
function [] = runExperiment( datasetName, imageExtension )
    %% ========== Step 0: Set program options and run initializations ==========
    %% Step 0.0: Get program options and parameters.
    options = SetParameters(datasetName);
    datasetFolder = [options.currentFolder '/input/' datasetName '/vocab/'];
    
    %% Step 0.1: Add relevant folders to path.
    addpath([options.currentFolder '/utilities']);
    addpath([options.currentFolder '/graphTools']);
    addpath([options.currentFolder '/vocabLearning']);
    
    %% Step 0.2: Specify name of the graph files and create output folders.
    graphFileName = [options.currentFolder '/graphs/' datasetName '.g'];
    resultFileName = [options.outputFolder '/' datasetName '.txt'];
    fp = fopen(graphFileName, 'w');
    
    % Create folder structures.
    if ~exist(options.processedFolder,'dir')
       mkdir(options.processedFolder); 
    end
    if ~exist(options.preDefinedFolder,'dir')
       mkdir(options.preDefinedFolder); 
    end
    if ~exist(options.testGraphFolder, 'dir')
       mkdir(options.testGraphFolder);
    end
    if ~exist(options.testOutputFolder, 'dir')
       mkdir(options.testOutputFolder); 
    end
    
    %% Step 0.3: Create initial data structures.
    fileNames = fuf([datasetFolder '*', imageExtension], 1, 'detail');
    trainingFileNames = fileNames;
    nodeCounter = 0;
    allNodes = cell(options.maxNumberOfFeatures,3);
    
    %% ========== Step 1: Pre-process the data (extract first level nodes, surpress weak responses) ==========
    if options.learnVocabulary
        %% Step 1.0: Downsample the image if it is too big.
        for fileItr = 1:size(trainingFileNames,1) 
            img = imread(trainingFileNames{fileItr});
            [~, fileName, ~] = fileparts(trainingFileNames{fileItr});
            if max(size(img)) > options.maxImageDim
               img = imresize(img, options.maxImageDim/max(size(img)), 'bilinear'); 
            end
            imwrite(img, [options.processedFolder '/' fileName '.png']);
        end
        
        %% Step 1.1: Extract a set of features from the input image.
        for fileItr = 1:size(trainingFileNames,1) 
            nodes = getNodes(img, options);
            % Keep nodes in the array.
            allNodes((nodeCounter + 1):(nodeCounter + size(nodes,1)), 1:2) = nodes;
            % Assign nodes their image ids.
            imageIds = ones(size(nodes,1), 1)*fileItr;
            allNodes((nodeCounter + 1):(nodeCounter + size(nodes,1)), 3) = ...
                                            mat2cell(imageIds, ones(size(imageIds)));
            % Increment node counter.
            nodeCounter = nodeCounter + size(nodes,1);
        end
        % Trim allNodes array to get rid of empty rows.
        allNodes = allNodes(1:nodeCounter,:);
    end
    %% ========== Step 2: Create first-level object graphs, and print them to a file. ==========
    % Upper level graphs are created in the main loop, as explained in
    % paper.
    if options.learnVocabulary
        %% Step 2.1: Get first-level object graph edges.
        [modes, edges, leafNodeAdjArr] = extractEdges(allNodes, [], [], options, 1, datasetName, []);

        %% Step 2.2: Print the object graphs to a file.
        imageIds = cell2mat(allNodes(:,3));
        numberOfImages = max(imageIds);
        nodeOffset = 0;
        for imageItr = 1:numberOfImages
            % Get only nodes and edges belonging to this image.
            imageNodeIdx = find(imageIds==imageItr);
            firstNodesOfEdges = edges(:,1);
            imageEdgeIdx = ismember(firstNodesOfEdges, imageNodeIdx);
            imageNodes = allNodes(imageIds==imageItr,:);
            imageEdges = edges(imageEdgeIdx, :);
            imageEdges(:,1:2) = imageEdges(:,1:2) - nodeOffset;
            
            % Print graph belonging to this image.
            fprintf(fp, 'XP\n');
            printGraphToFile(fp, imageNodes(:,1), imageEdges, true);
            nodeOffset = nodeOffset + size(imageNodes,1);
        end
        fclose(fp);
    end
    %% ========== Step 3: Create compositional vocabulary (Main loop in algorithm 1 of paper). ==========
    if options.learnVocabulary
        [vocabulary, mainGraph, modes] = learnVocabulary(allNodes, edges, modes, leafNodeAdjArr, graphFileName, ...
                                        resultFileName, options, trainingFileNames, datasetName);
        save([options.currentFolder '/output/' datasetName '_vb.mat'], 'vocabulary', 'mainGraph', 'modes', 'leafNodeAdjArr');
    else
        load([options.currentFolder '/output/' datasetName '_vb.mat'], 'vocabulary', 'mainGraph', 'modes', 'leafNodeAdjArr');
    end
    %% ========== Step 4: Run inference for all test images with the learned vocabulary. ==========
    if options.testImages
        testFileNames = fileNames;
        %% Step 4.1: Create files for pre-defined substructures ( compositions from voc. at each level)
        preparePreDefinedFiles(options.preDefinedFolder, vocabulary);
        
        %% Step 4.2: Run inference on each test image.
        for testImgItr = 1:size(testFileNames,1)
            singleTestImage(testFileNames{testImgItr}, options, options.currentFolder);
        end
        
        %% Step 4.3: Run image retrieval tests.
%        runRetrievalTests(options.currentFolder, options.datasetName, options.testGraphFolder, 3);
    end
end

