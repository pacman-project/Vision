%> Name: runVocabularyLearning
%>
%> Description: The entry function to CHOP code. Calls CHOP with
%> desired parameters on specified dataset. For each image in the dataset,
%> we create a graph out of level 1 gabor filter responses, and form
%> edges. The overall graph that consists of n (image count)
%> non-overlapping graphs is called "main graph" (how creative). The main
%> graph is compressed to obtain parts for the next level. Overlapping part
%> realizations in this new level are processed (so low-rank parts are
%> inhibited). The remaining part realizations contribute to the main graph 
%> of the next level. This procedure lasts until no new parts are found.
%>
%> @param datasetName Name of the dataset to work on. 
%> @param imageExtension The extension of the files to work on. Examples
%> include '.jpg', '.png', '_crop.png'...
%> @param imageExtension The ground truth image extension for each file.
%> For example, for a training image "swan.png", the GT image can have the 
%> name "swan_gt.png", where imageExtension is '.png', and gtImageExtension 
%> is '_gt.png'.
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 18.11.2013
%> Ver 1.1 on 05.12.2013 Various parameter additions, 'mode' changes
%> Ver 1.2 on 12.01.2014 Comment changes for unified code look
%> Ver 1.3 on 12.01.2014 Timing is added by Mete
%> Ver 1.4 on 17.02.2014 GT processing added.
%> Ver 1.5 on 16.06.2014 New comments, simplification on code.
function [] = runVocabularyLearning( datasetName, imageExtension, gtImageExtension )
    %% ========== Step 0: Set program options and run initializations ==========
    %% Step 0.0: Get program options and parameters.
    options = SetParameters(datasetName, true);
    datasetFolder = [options.currentFolder '/input/' datasetName '/vocab/'];
    gtFolder = [options.currentFolder '/input/' datasetName '/gt/'];
    processedFolder = options.processedFolder;
    processedGTFolder = options.processedGTFolder;
    
    % Open threads for parallel processing.
    if options.parallelProcessing && usejava('jvm')
        if options.parpoolFlag
            p = gcp('nocreate');
            if ~isempty(p)
                delete(p);
            end
            parpool(options.numberOfThreads);
        else
            s = matlabpool('size');
            if s>0
               matlabpool close; 
            end
            matlabpool('open', options.numberOfThreads);
        end
    end
    
    if options.learnVocabulary
        % Create the folder structure required.
        createFolders(options);

        %% Step 0.1: Create initial data structures.
        try
            fileNames = fuf([datasetFolder '*', imageExtension], 1, 'detail');
        catch
            display(['Unable to find training images. Possibly you forgot to put the images under ./input/' datasetName '/vocab/ or defined the extension wrong.']); 
            return;
        end
        trainingFileNames = fileNames;
        if isempty(trainingFileNames)
           display(['Unable to find training images. Possibly you forgot to put the images under ./input/' datasetName '/vocab/ or defined the extension wrong.']);  
        end

        %% Step 0.2: Allocate space to keep names of corresponding gt files.
        gtFileNames = cell(numel(trainingFileNames),1);
        
        %% ========== Step 1: Pre-process the data (extract first level nodes, surpress weak responses) ==========
        %% Step 1.0: Pre-processing. Downsample images, and associate them with their ground truth.
        maxImageDim = options.maxImageDim;
        for fileItr = 1:size(trainingFileNames,1) 
            % Read image and downsample it.
            img = imread(trainingFileNames{fileItr});
            [~, fileName, ~] = fileparts(trainingFileNames{fileItr});
            if max(size(img)) > maxImageDim
               img = imresize(img, maxImageDim/max(size(img)), 'bilinear'); 
            end
            
            % If img is binary, we can save it as a binary png.
            if size(img,3) == 1 && max(max(max(img))) == 1
                img = img > 0;
            end
            
            % Save image into processed folder.
            fileNameNew = fileName;
            if exist([processedFolder '/' fileName '.png'], 'file')
                fileId = 1;
                fileNameNew = [fileName '_' num2str(fileId)];
                while exist([processedFolder '/' fileNameNew '.png'], 'file')
                    fileId = fileId + 1;
                    fileNameNew = [fileName '_' num2str(fileId)];
                end
            end
            imwrite(img, [processedFolder '/' fileNameNew '.png']);
            
            % Switch file names with those copied.
            trainingFileNames(fileItr) = {[processedFolder '/' fileNameNew '.png']};

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
        display('..... Level 1 Node Extraction started. This may take a while.');
        allNodes = cell(size(trainingFileNames,1),1);
        allNodeActivations = cell(size(trainingFileNames,1),1);
        smoothedFolder = options.smoothedFolder;
        
        % If the images have been processed earlier, no need for us to do
        % this again.
        if exist([options.currentFolder '/outputNodes/' options.datasetName '/nodes.mat'], 'file')
           load([options.currentFolder '/outputNodes/' options.datasetName '/nodes.mat']);
        else
            parfor fileItr = 1:size(trainingFileNames,1)
                [~, fileName, ~] = fileparts(trainingFileNames{fileItr});
                img = imread([processedFolder '/' fileName '.png']);

    %            edgeImg = edge(img, 'canny', [0.05 0.1], 1);
                % Get the Level 1 features.
                [nodes, smoothedImg, nodeActivations, smoothActivationImg, responseImgs, edgeImg] = getNodes(img, gtFileNames{fileItr}, options);

                % Generate an imperfect backprojection from obtained peaks.
    %             level1Nodes = cell2mat(nodes);
     %            level1Nodes = [level1Nodes(:,[1,4,5]), double(nodeActivations)];
    %             img1 = obtainPoE(level1Nodes, size(edgeImg), options);
    %             
    %             % Generate a perfect backprojection.
    %             [activationImg, nodeIdImg] = max(responseImgs, [], 3);
    %             activationImg(edgeImg == 0) = 0;
    %             nodeIdImg(edgeImg == 0) = 0;
    %             peaks = find(nodeIdImg);
    %             [x,y] = ind2sub(size(nodeIdImg), peaks);
    %             level1Nodes = [nodeIdImg(peaks), x, y, activationImg(peaks)];
    %             img2 = obtainPoE(level1Nodes, size(nodeIdImg), options);

                % Keep nodes in the array.
                allNodes(fileItr) = {nodes};
                allNodeActivations(fileItr) = {nodeActivations};

                % Save smoothed image.
                imwrite(smoothedImg, [smoothedFolder '/' fileName '_peaks.png']);
                imwrite(edgeImg, [smoothedFolder '/' fileName '_edgeImg.png']);
                imwrite(smoothActivationImg, [smoothedFolder '/' fileName '.png']);
    %  %           imwrite(img1, [smoothedFolder '/' fileName '_peakPoE.png']);
    % %            imwrite(img2, [smoothedFolder '/' fileName '_perfectPoE.png']);
                parsave([smoothedFolder '/' fileName '_responseImgs.mat'], responseImgs, 'responseImgs');
            end
            if ~exist([options.currentFolder '/outputNodes/' options.datasetName], 'dir')
                mkdir([options.currentFolder '/outputNodes/' options.datasetName]);
            end
            save([options.currentFolder '/outputNodes/' options.datasetName '/nodes.mat'], 'allNodes', 'allNodeActivations', '-v7');
        end
        
        % Calculate image size.
        [~, fileName, ~] = fileparts(trainingFileNames{1});
        img = imread([processedFolder '/' fileName '.png']);
        imageSize = size(img);
        options.imageSize = imageSize(1:2);

        % Reorder images based on their node count. This helps in
        % efficient parallelization. 
        nodeCounts = cellfun(@(x) size(x,1), allNodes);
        [~, sortedImageIdx] = sort(nodeCounts, 'descend');
        trainingFileNames = trainingFileNames(sortedImageIdx);
        allNodes = allNodes(sortedImageIdx);
        allNodeActivations = allNodeActivations(sortedImageIdx);

        %% Learn category/pose of each image and and correct their order 
        % based on the sort order of the images.
        categoryArr = cell(numel(trainingFileNames),1);
        poseArr = zeros(numel(trainingFileNames),1);
        for fileItr = 1:numel(fileNames)
            % Get category label.
            fullName = fileNames{fileItr};
            [~, fileName, ext] = fileparts(fileNames{fileItr});
            strLength = numel([options.currentFolder '/input/' options.datasetName '/vocab/']);
            fileNameLength = numel(ext) + numel(fileName) + 1; % 1 for '/' (folder delimiter)
            if numel(fullName) >= strLength + fileNameLength
                categoryStr = fullName(1, (strLength+1):(end-fileNameLength));
                categoryStrSepIdx = strfind(categoryStr, '/');
                if ~isempty(categoryStrSepIdx)
                    categoryStr = categoryStr(:, 1:(categoryStrSepIdx(1)-1));
                end
                categoryArr(fileItr) = {categoryStr};
            else
                categoryArr(fileItr) = {''};
            end

            % Get pose info.
            poseStartIdx = strfind(fileName, '_r');
            if ~isempty(poseStartIdx)
                poseIdx = poseStartIdx + numel('_r');
                poseArr(fileItr) = sscanf(fileName(1, poseIdx:end), '%d');
            end
        end
        categoryArr = categoryArr(sortedImageIdx);
        poseArr = poseArr(sortedImageIdx); %#ok<NASGU>
       
        % Set image ids for each node.
        imageIds = cell(size(allNodes,1),1);
        for fileItr = 1:size(allNodes,1)
            imageIds(fileItr) = {num2cell(repmat(fileItr, size(allNodes{fileItr},1), 1))};
        end
        
        % Convert all node data into a single matrix.
        allNodes = cat(1, allNodes{:});
        allNodeActivations = cat(1, allNodeActivations{:});
        imageIds = cat(1, imageIds{:});
        leafNodes = [allNodes, imageIds];
        leafNodes = int32(cell2mat(leafNodes));
        leafNodeCoords = leafNodes(:,4:5);
        leafNodes = leafNodes(:, [1:3 6]);
        clear allNodes; 
        clear fileNames;
        
        % Learn node signs based on whether or not they are from the
        % background class. The ones which are from the background class
        % will have negative signs, while all others have positive signs.
        % The effect is that in Subdue, subgraphs that compress positive
        % graphs but not negative ones will be favoured.
        imageSigns = ~strcmp(categoryArr, options.backgroundClass);
        leafNodeSigns = imageSigns(cell2mat(imageIds));
        clear imageIds;
        
        % Transform category array into an index-based one (having numbers
        % instead of category strings). The category string labels is saved
        % in categoryNames.
        [~, categoryNames, categoryArrIdx] = unique(categoryArr, 'stable'); 
        categoryNames = categoryArr(categoryNames); 
        options.numberOfCategories = numel(categoryNames);
        clear categoryArr;
        
        %% ========== Step 2: Create first-level object graphs, and print them to a file. ==========
        [vocabLevel, vocabLevelDistributions, graphLevel] = generateLevels(leafNodes, leafNodeCoords, allNodeActivations, leafNodeSigns, options);
        clear allNodeActivations;

        %% Learn edge-based distance matrix once and for all.
        [edgeIdMatrix, edgeDistanceMatrix, edgeCoords, edgeLogMin, edgeLogRange ] = findEdgeDistanceMatrix(options.receptiveFieldSize, options.distType, options.edgeSimilarityAllowed);
        options.edgeIdMatrix = edgeIdMatrix;
        options.edgeDistanceMatrix = edgeDistanceMatrix;
        options.edgeCoords = edgeCoords;
        options.edgeLogMin = edgeLogMin;
        options.edgeLogRange = edgeLogRange;
        clear edgeIdMatrix edgeDistanceMatrix edgeCoords edgeLogMin edgeLogRange;
        
        %% Step 2.1: Get first-level object graph edges.
        [graphLevel] = extractEdges(graphLevel, [], options, 1);
        
        %% Here, we bring back statistical learning with mean/variance.
        [modes, modeProbArr] = learnModes(graphLevel, options.edgeCoords, options.edgeIdMatrix, options.datasetName, 1, options.currentFolder, ~options.fastStatLearning && options.debug);
        graphLevel = assignEdgeLabels(graphLevel, modes, modeProbArr, options.edgeCoords, 1, options.debugFolder);
        
        %% ========== Step 3: Create compositional vocabulary (Main loop in algorithm 1 of ECCV 2014 paper). ==========
        tr_s_time=tic;  
        save([options.currentFolder '/output/' datasetName '/export.mat'], 'trainingFileNames', 'categoryNames', 'categoryArrIdx', 'poseArr', '-v7');
        save([options.currentFolder '/output/' datasetName '/vb.mat'], 'trainingFileNames', 'categoryNames', '-v7');
        
        %% Learn vocabulary!
        learnVocabulary(vocabLevel, vocabLevelDistributions, graphLevel, leafNodes, leafNodeCoords, options, trainingFileNames, modes, modeProbArr);
        
        % Stop counting time.
        tr_stop_time=toc(tr_s_time); %#ok<NASGU>
    end
    
    % Close thread pool if opened.
    if options.parallelProcessing
        if ~options.parpoolFlag
            matlabpool close;
        end
    end
end

function parsave(fname, data, dataName)
     save(fname)
end