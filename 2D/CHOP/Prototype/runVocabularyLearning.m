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
function [] = runVocabularyLearning( datasetName, imageExtension )
    %% ========== Step 0: Set program options and run initializations ==========
    %% Step 0.0: Get program options and parameters.
    options = SetParameters(datasetName);
    datasetFolder = [options.currentFolder '/input/' datasetName '/vocab/'];
    
    if options.learnVocabulary
        %% Step 0.2: Create initial data structures.
        fileNames = fuf([datasetFolder '*', imageExtension], 1, 'detail');
        trainingFileNames = fileNames;
        nodeCounter = 0;
        allNodes = cell(options.maxNumberOfFeatures,6);
    
        %% ========== Step 1: Pre-process the data (extract first level nodes, surpress weak responses) ==========
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
            [~, fileName, ~] = fileparts(trainingFileNames{fileItr});
            img = imread([options.processedFolder '/' fileName '.png']);
            nodes = getNodes(img, options);
            % Keep nodes in the array.
            allNodes((nodeCounter + 1):(nodeCounter + size(nodes,1)), 1:2) = nodes;
            % Assign nodes their image ids.
            imageIds = ones(size(nodes,1), 1)*fileItr;
            allNodes((nodeCounter + 1):(nodeCounter + size(nodes,1)), 3) = ...
                                            mat2cell(imageIds, ones(size(imageIds)));
            allNodes((nodeCounter + 1):(nodeCounter + size(nodes,1)), 4) = ...
                                            mat2cell(zeros(size(imageIds)), ones(size(imageIds)));
            allNodes((nodeCounter + 1):(nodeCounter + size(nodes,1)), 5) = ...
                                            mat2cell(ones(size(imageIds)), ones(size(imageIds)));
            % Increment node counter.
            nodeCounter = nodeCounter + size(nodes,1);
        end
        % Trim allNodes array to get rid of empty rows.
        allNodes = allNodes(1:nodeCounter,:);
        
        %% Step 1.2: If receptive field is used, nodes will be repeated.
        % (so that each node set corresponds to a different receptive
        % field)
        leafNodes = allNodes;
        [allNodes, ~] = getReceptiveFieldNodes(allNodes, 1, options);
        
        %% ========== Step 2: Create first-level object graphs, and print them to a file. ==========
        [vocabLevel, graphLevel] = generateLevels(allNodes, leafNodes, options);
        
        %% Step 2.1: Get first-level object graph edges.
        mainGraph = {graphLevel};
        [modes, highLevelModes, mainGraph, leafNodeAdjArr] = extractEdges(mainGraph, [], options, 1, datasetName, [], []);
        graphLevel = mainGraph{1};
        
        %% ========== Step 3: Create compositional vocabulary (Main loop in algorithm 1 of paper). ==========
        tr_s_time=tic;    

        [vocabulary, mainGraph, modes, highLevelModes] = learnVocabulary(vocabLevel, graphLevel, leafNodes(:,1:3), modes, highLevelModes, leafNodeAdjArr, ...
                                        options, trainingFileNames, datasetName);
                                    
        tr_stop_time=toc(tr_s_time)
 

        save([options.currentFolder '/output/' datasetName '/' datasetName '_trtime.mat'], 'tr_stop_time');

        save([options.currentFolder '/output/' datasetName '/' datasetName '_vb.mat'], 'vocabulary', 'mainGraph', 'modes', 'highLevelModes', 'leafNodeAdjArr', 'leafNodes', 'fileNames');
    end
end

