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
    if options.parallelProcessing
        s = matlabpool('size');
        if s>0
           matlabpool close; 
        end
        matlabpool('open', options.numberOfThreads);
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
            imwrite(img, [processedFolder '/' fileName '.png']);

            % Switch file names with those copied.
            trainingFileNames(fileItr) = {[processedFolder '/' fileName '.png']};

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
        smoothedFolder = options.smoothedFolder;
        parfor fileItr = 1:size(trainingFileNames,1)
            [~, fileName, ~] = fileparts(trainingFileNames{fileItr});
            img = imread([processedFolder '/' fileName '.png']);
            
            % Get the Level 1 features.
            [nodes, smoothedImg] = getNodes(img, gtFileNames{fileItr}, options);

            % Keep nodes in the array.
            allNodes(fileItr) = {nodes};

            % Save smoothed image.
            imwrite(smoothedImg, [smoothedFolder '/' fileName '.png']);
        end

        % Reorder images based on their node count. This helps in
        % efficient parallelization. 
        nodeCounts = cellfun(@(x) size(x,1), allNodes);
        [~, sortedImageIdx] = sort(nodeCounts, 'descend');
        trainingFileNames = trainingFileNames(sortedImageIdx);
        allNodes = allNodes(sortedImageIdx);

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
        imageIds = cat(1, imageIds{:});
        leafNodes = [allNodes, imageIds];
        
        % Learn node signs based on whether or not they are from the
        % background class. The ones which are from the background class
        % will have negative signs, while all others have positive signs.
        % The effect is that in Subdue, subgraphs that compress positive
        % graphs but not negative ones will be favoured.
        imageSigns = ~strcmp(categoryArr, options.backgroundClass);
        leafNodeSigns = imageSigns(cell2mat(imageIds));
        
        %% ========== Step 2: Create first-level object graphs, and print them to a file. ==========
        [vocabLevel, graphLevel] = generateLevels(leafNodes, leafNodeSigns, options);
        
        %% Step 2.1: Get first-level object graph edges.
        mainGraph = {graphLevel};
        [modes, mainGraph] = extractEdges(mainGraph, options, 1, []);
        graphLevel = mainGraph{1};
        
        %% ========== Step 3: Create compositional vocabulary (Main loop in algorithm 1 of ECCV 2014 paper). ==========
        tr_s_time=tic;  
        [vocabulary, redundantVocabulary, mainGraph, modes, distanceMatrices] = learnVocabulary(vocabLevel, graphLevel, leafNodes(:,1:3), modes, ...
                                        options, trainingFileNames); %#ok<NASGU,ASGLU>
        tr_stop_time=toc(tr_s_time); %#ok<NASGU>
        
        % Export realizations into easily-readable arrays.
        exportArr = exportRealizations(mainGraph); %#ok<NASGU>
        
        % Transform category array into an index-based one (having numbers
        % instead of category strings). The category string labels is saved
        % in categoryNames.
        [~, categoryNames, categoryArrIdx] = unique(categoryArr, 'stable'); %#ok<NASGU>
        categoryNames = categoryArr(categoryNames); %#ok<NASGU>
        
        % Print everything to files.
        save([options.currentFolder '/output/' datasetName '/trtime.mat'], 'tr_stop_time');
        save([options.currentFolder '/output/' datasetName '/vb.mat'], 'vocabulary', 'redundantVocabulary', 'modes', 'trainingFileNames', 'categoryNames');
        % categoryArr is kept for backward-compatibility. It will be
        % removed in further releases.
        save([options.currentFolder '/output/' datasetName '/export.mat'], 'trainingFileNames', 'exportArr', 'categoryArr', 'categoryArrIdx', 'poseArr'); 
    end
    
    % Close thread pool if opened.
    if options.parallelProcessing
        matlabpool close;
    end
end