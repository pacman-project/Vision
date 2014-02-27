%> Name: runVocabularyLearning
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
%> Ver 1.3 on 12.01.2014 Timing is added by Mete
%> Ver 1.4 on 17.02.2014 GT processing added.
function [] = runVocabularyLearning( datasetName, imageExtension, gtImageExtension )
    %% ========== Step 0: Set program options and run initializations ==========
    %% Step 0.0: Get program options and parameters.
    options = SetParameters(datasetName);
    datasetFolder = [options.currentFolder '/input/' datasetName '/vocab/'];
    gtFolder = [options.currentFolder '/input/' datasetName '/gt/'];
    processedFolder = options.processedFolder;
    processedGTFolder = options.processedGTFolder;
    
    % Open threads for parallel processing.
    matlabpool('open', options.numberOfThreads);
    
    if options.learnVocabulary
        %% Step 0.1: Create initial data structures.
        fileNames = fuf([datasetFolder '*', imageExtension], 1, 'detail');
        trainingFileNames = fileNames;
        
        %% Step 0.2: Allocate space to keep names of corresponding gt files.
        gtFileNames = cell(numel(trainingFileNames),1);
        %% ========== Step 1: Pre-process the data (extract first level nodes, surpress weak responses) ==========
        %% Step 1.0: Downsample the image if it is too big.
        maxImageDim = options.maxImageDim;
        parfor fileItr = 1:size(trainingFileNames,1) 
            % Read image and downsample it.
            img = imread(trainingFileNames{fileItr});
            [~, fileName, ~] = fileparts(trainingFileNames{fileItr});
            if max(size(img)) > maxImageDim
               img = imresize(img, maxImageDim/max(size(img)), 'bilinear'); 
            end
            imwrite(img, [processedFolder '/' fileName '.png']);
            
            % If gt file exists, write it to the processed gt folder. In
            % addition, we keep track of the name of gt file corresponding
            % to the image read.
            gtFile = [gtFolder fileName gtImageExtension];
            if exist(gtFile, 'file')
                gtImg = imread(gtFile);
                if size(gtImg,3)>2
                    gtImg = rgb2gray(gtImg(:,:,1:3))>0;
                else
                    gtImg = gtImg(:,:,1)>0;
                end
                gtImg = imresize(gtImg, [size(img,1), size(img,2)], 'bilinear');
                gtFileNames(fileItr) = {[processedGTFolder '/' fileName '.png']};
                imwrite(gtImg, [processedGTFolder '/' fileName '.png']);
            end
        end
        
        %% Step 1.1: Extract a set of features from the input images.
        allNodes = cell(size(trainingFileNames,1),1);
        parfor fileItr = 1:size(trainingFileNames,1)
            [~, fileName, ~] = fileparts(trainingFileNames{fileItr});
            img = imread([processedFolder '/' fileName '.png']);
            nodes = getNodes(img, gtFileNames{fileItr}, options);
            
            % Keep nodes in the array.
            allNodes(fileItr) = {nodes};
        end
        imageIds = cell(size(allNodes,1),1);
        for fileItr = 1:size(allNodes,1)
            imageIds(fileItr) = {num2cell(repmat(fileItr, size(allNodes{fileItr},1), 1))};
        end
        
        allNodes = cat(1, allNodes{:});
        imageIds = cat(1, imageIds{:});
        allNodes = [allNodes, imageIds];
        
        %% Step 1.2: If receptive field is used, nodes will be repeated.
        % (so that each node set corresponds to a different receptive
        % field)
        leafNodes = allNodes;
        
        %% ========== Step 2: Create first-level object graphs, and print them to a file. ==========
        [vocabLevel, graphLevel] = generateLevels(allNodes, options);
        
        %% Step 2.1: Get first-level object graph edges.
        mainGraph = {graphLevel};
        [modes, highLevelModes, mainGraph] = extractEdges(mainGraph, options, 1, [], []);
        graphLevel = mainGraph{1};
        
        %% ========== Step 3: Create compositional vocabulary (Main loop in algorithm 1 of paper). ==========
        tr_s_time=tic;  
        [vocabulary, mainGraph, modes, highLevelModes] = learnVocabulary(vocabLevel, graphLevel, leafNodes(:,1:3), modes, highLevelModes, ...
                                        options, trainingFileNames, datasetName);
        tr_stop_time=toc(tr_s_time);
        save([options.currentFolder '/output/' datasetName '/' datasetName '_trtime.mat'], 'tr_stop_time');
        save([options.currentFolder '/output/' datasetName '/' datasetName '_vb.mat'], 'vocabulary', 'mainGraph', 'modes', 'highLevelModes', 'leafNodes', 'fileNames');
    end
    matlabpool close;
end

