%> Name: learnVocabulary
%>
%> Description: Given the graph of first-level responses and their
%> relations, this function extracts a hierarchical vocabulary by detecting
%> structural patterns inherent in the input graph. The pattern discovery
%> phase is hierarchic, and continues until a the graph is compressed into
%> a single node.
%>
%> @param vocabLevel The first level of vocabulary, level 1 parts.
%> @param graphLevel Object graphs level 1.
%> @param options Program options
%> @param fileList input image name list.
%>
%> @retval vocabulary The hierarchic vocabulary learnt from the data.
%> @retval mainGraph The hierarchic object graphs.
%> @retval distanceMatrices Distance matrix of compositions.
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 26.11.2013
%> Ver 1.1 on 02.12.2013 Completed unlimited number of graph generation.
%> Ver 1.2 on 12.01.2014 Commentary changes, for unified code look.
%> Ver 1.3 on 28.01.2014 Mode calculation put in a separate file.
%> Ver 1.4 on 03.02.2014 Refactoring
   % Reduce memory consumption by writing all stuff to files,
   % clearing all, and then returning back to computation.
   addpath(genpath([pwd '/utilities']));
   if exist([pwd '/Workspace.mat'], 'file')
        disp('Loading workspace from the previous layer.');
        load([pwd '/Workspace.mat']);
   end
   
  % Open/close matlabpool to save memory. 
   if options.parallelProcessing
       s = matlabpool('size');
       if s>0
          matlabpool close; 
       end
       matlabpool('open', options.numberOfThreads);
   end
   % Increase level itr.
   levelItr = levelItr + 1; %#ok<*NASGU>
   
   % If this is a fresh Matlab, modify paths.
   % ========== PATH FOLDER ADDITION ==========
   if options.restartFlag
       w = warning('off', 'all');
       addpath(genpath([options.currentFolder '/demo']));
       addpath(genpath([options.currentFolder '/graphTools']));
       addpath(genpath([options.currentFolder '/vocabLearning']));
       addpath(genpath([options.currentFolder '/inference']));
       addpath(genpath([options.currentFolder '/categorization']));
       warning(w);
   end
   
   % Try to restore break points.
   try %#ok<TRYNC>
   dbstop(bInfo);
   end
   
   % Obtain the pre-set threshold for this level, if there is one.
   if levelItr > (numel(options.subdue.presetThresholds) + 1)
       presetThreshold = options.subdue.threshold;
       options.optimizationFlag = orgOptimizationFlag;
   else
       presetThreshold = options.subdue.presetThresholds(levelItr-1);
       options.optimizationFlag = false;
   end

   %% Step 2.1: Run knowledge discovery to learn frequent compositions.
   [vocabLevel, graphLevel, ~, isSupervisedSelectionRunning, previousAccuracy] = discoverSubs(vocabLevel, graphLevel, newDistanceMatrix,...
       options, presetThreshold, levelItr-1, supervisedSelectionFlag, isSupervisedSelectionRunning, previousAccuracy, level1Coords); %#ok<*ASGLU,*NODEF>
   java.lang.System.gc();
   % Test if the mapping is correct (for debugging).
%        testSuccess = testMapping(vocabLevel, graphLevel, newDistanceMatrix, mainGraph{levelItr-1});

   % Open/close matlabpool to save memory.
   matlabpool close;
   matlabpool('open', options.numberOfThreads);

   %% If no new subs have been found, finish processing.
   if isempty(vocabLevel)
      % Write previous level's appearances to the output folder.
      vocabulary = vocabulary(1:(levelItr-1),:);
      vocabularyDistributions = vocabularyDistributions(1:(levelItr-1),:);
      mainGraph = mainGraph(1:(levelItr-1),:);
      allModes = allModes(1:(levelItr-1), :);
      distanceMatrices = distanceMatrices(1:(levelItr-1),:);
      modeProbs = modeProbs(1:(levelItr-1),:);
      save([options.currentFolder '/output/' options.datasetName '/mainGraph.mat'], 'mainGraph', '-v7.3'); 
      return; 
   end

   %% Assign realizations R of next graph level (l+1), and fill in their bookkeeping info.
   previousLevel = mainGraph{levelItr-1};
   graphLevel = fillBasicInfo(previousLevel, graphLevel, levelItr, options);
   java.lang.System.gc();

   %% Calculate statistics from this graph.
   display('........ Before we apply inhibition, estimating statistics..');
   [avgShareability, avgCoverage] = saveStats(vocabLevel, graphLevel, leafNodeCoords, maxCoverageVals, numberOfImages, options, 'preInhibition', levelItr);
   display(['........ Average Coverage: ' num2str(avgCoverage) ', average shareability of compositions: ' num2str(avgShareability) ' percent.']); 

   %% If category level is reached, we reduce the number of desired nodes substantially. 
   nodeCoverages = zeros(numel(graphLevel),1);
   imageIds = [graphLevel.imageId];
   remainingImages = unique([graphLevel.imageId]);
   imageCoverages = zeros(max(imageIds), 1);
   stopFlag = false;
   for nodeItr = 1:numel(graphLevel)
        nodeCoverages(nodeItr) = (numel(graphLevel(nodeItr).leafNodes) / nnz(level1Coords(:,1) == graphLevel(nodeItr).imageId)) / avgCoverage;
        imageCoverages(graphLevel(nodeItr).imageId) = max(imageCoverages(graphLevel(nodeItr).imageId), nodeCoverages(nodeItr));
   end
   if levelItr == options.categoryLevel || mean(imageCoverages(remainingImages)) >= options.categoryLevelCoverage
        stopFlag = true;
        options.categoryLevel = levelItr;
        options.stopFlag = true;
   end

   %% In order to do proper visualization, we learn precise positionings of children for every vocabulary node.
   display('........ Learning sub-part label and position distributions.');
   [vocabLevel, vocabLevelDistributions] = learnChildDistributions(vocabLevel, graphLevel, mainGraph{levelItr-1}, levelItr, options);
   java.lang.System.gc();

   %% Calculate activations for every part realization.
   display('........ Calculating activations..');
   graphLevel = calculateActivations(vocabLevel, vocabulary, graphLevel, mainGraph, level1Coords, options, levelItr);
   java.lang.System.gc();

   %% Remove low-coverage nodes when we're at category layer.
   if stopFlag
        graphLevel = graphLevel(nodeCoverages >= min(options.categoryLevelCoverage, imageCoverages(imageIds)));
   end
   
   %% Inhibition! We process the current level to eliminate some of the nodes in the final graph.
   % The rules here are explained in the paper. Basically, each node
   % should introduce a novelty (cover a new area of the image). If it
   % fails to do so, the one which has a lower mdl ranking is
   % discarded. *Natural selection*
   if options.noveltyThr > 0
       display('........ Applying inhibition..');
       graphLevel=applyTestInhibition(graphLevel, options);
       display(['........ Inhibition applied with novelty thr: ' num2str(options.noveltyThr) '.']);
       display(['........ We have ' num2str(numel(graphLevel)) ' realizations of ' num2str(numel(unique([graphLevel.labelId]))) ' compositions.']);
   else
       display('........ No inhibition applied..');
   end

   % Open/close matlabpool to save memory.
   matlabpool close;
   matlabpool('open', options.numberOfThreads);

   %% Post-process graphLevel, vocabularyLevel to remove non-existent parts from vocabLevel.
   % In addition, we re-assign the node ids in graphLevel.
  display('........ Calculating distance matrix among the vocabulary nodes (in parallel)..');
   vocabularyDistributions{levelItr} = vocabLevelDistributions;
  [vocabLevel, graphLevel, eucDistanceMatrix] = postProcessParts(vocabLevel, graphLevel, vocabulary, vocabularyDistributions, levelItr, options);       
  
  % Visualize learned OR Nodes.
  visualizeORNodes( vocabLevel, levelItr, options);
  
  % Update distance matrices.
  newDistanceMatrix = inf(size(eucDistanceMatrix,1), size(eucDistanceMatrix,1), 'single');
  newDistanceMatrix(1:(size(eucDistanceMatrix,1)+1):size(eucDistanceMatrix,1)*size(eucDistanceMatrix,1)) = 0;
  distanceMatrices{levelItr} = newDistanceMatrix;
   java.lang.System.gc();

   %% We process graphLevel's labelIds to reflect updated labels (OR Node Labels).
   vocabNodeLabels = [vocabLevel.label];
   updatedLabelIds = num2cell(vocabNodeLabels([graphLevel.labelId]));
   [graphLevel.labelId] = deal(updatedLabelIds{:});
   mainGraph{levelItr} = graphLevel;

   %% As an initial stage of inhibition, we downsample the responses, 
   % and apply max pooling. Please note that turning this off also means no
   % reduction in resolution for higher layers, therefore non-growing
   % receptive fields.
   display(['........ Applying pooling on ' num2str(numel(graphLevel)) ' realizations belonging to ' num2str(max([vocabLevel.label])) ' compositions.']);
   graphLevel = applyPooling(graphLevel, options.poolDim, options.poolFlag, ~ismember(levelItr, options.noPoolingLayers));
   java.lang.System.gc();

   %% Calculate statistics from this graph.
   display('........ Estimating post-inhibition statistics..');
   [avgShareability, avgCoverage] = saveStats(vocabLevel, graphLevel, leafNodeCoords, maxCoverageVals, numberOfImages, options, 'postInhibition', levelItr);

   %% Experimenting. After some point, we need to convert to centroid-based edge creation, no matter what.
   if avgCoverage < options.minContinuityCoverage && options.edgeChangeLevel == -1 && ~strcmp(options.edgeType, 'centroid') || levelItr == options.maxEdgeChangeLevel
       options.edgeType = 'centroid';
       display('........ Switching to -centroid- type edges!');
       options.edgeChangeLevel = levelItr;
   end

   % display debugging info.
   display(['........ Remaining: ' num2str(numel(graphLevel)) ' realizations belonging to ' num2str(max([vocabLevel.label])) ' compositions.']);
   display(['........ Average Coverage: ' num2str(avgCoverage) ', average shareability of compositions: ' num2str(avgShareability) ' percent.']); 

   %% Step 2.5: Create the parent relationships between current level and previous level.
   vocabulary = mergeIntoGraph(vocabulary, vocabLevel, leafNodes, levelItr, 0);
   mainGraph = mergeIntoGraph(mainGraph, graphLevel, leafNodes, levelItr, 1);

      %% If we're at category level, our parts are supposed to cover most of the object. We can enforce this.
    if levelItr == options.categoryLevel - 1
         options.edgeNoveltyThr = options.categoryLevelEdgeNoveltyThr;
         options.maxNodeDegree = options.maxNodeDegree * 2; 
    end
   
   %% Step 2.6: Create object graphs G_(l+1) for the next level, l+1.
   % Extract the edges between new realizations to form the new object graphs.
   [mainGraph] = extractEdges(mainGraph, options, levelItr);
   graphLevel = mainGraph{levelItr};
   java.lang.System.gc();

   %% Here, we bring back statistical learning with mean/variance.
   [modes, modeProbArr] = learnModes(graphLevel, options.edgeCoords, options.edgeIdMatrix, options.datasetName, levelItr, options.currentFolder, ~options.fastStatLearning && options.debug);
   graphLevel = assignEdgeLabels(graphLevel, modes, modeProbArr, options.edgeCoords);
   mainGraph{levelItr} = graphLevel;
   allModes{levelItr} = modes;
   modeProbs{levelItr} = modeProbArr;
   java.lang.System.gc();

   %% Print distance matrices.
   if ~isempty(newDistanceMatrix)
      imwrite(eucDistanceMatrix, [options.currentFolder '/debug/' options.datasetName '/level' num2str(levelItr) '_dist.png']);
      imwrite(newDistanceMatrix, [options.currentFolder '/debug/' options.datasetName '/level' num2str(levelItr) '_orNodes.png']);
   end
   
   %% At this step, we perform imagination of the parts at this layer (for later use).
   vocabularyDistributions = projectVocabulary( vocabularyDistributions);

   %% If we've reached max number of layers, don't keep going forward.
   if levelItr == options.maxLevels || stopFlag
       vocabulary = vocabulary(1:(levelItr),:);
       vocabularyDistributions = vocabularyDistributions(1:(levelItr),:);
       allModes = allModes(1:(levelItr), :);
       mainGraph = mainGraph(1:(levelItr),:);
       distanceMatrices = distanceMatrices(1:(levelItr),:);
       modeProbs = modeProbs(1:(levelItr),:);
   end
   
   %% We're exporting output here. This helps us to perform 
   % Export realizations into easily-readable arrays.
   display('Writing output to files.');
   [exportArr, activationArr] = exportRealizations(mainGraph); %#ok<ASGLU>
   save([options.currentFolder '/output/' options.datasetName '/export.mat'], 'exportArr', 'activationArr', '-append', '-v7.3'); 
   clear exportArr  activationArr precisePositions;

   % Print everything to files.
   save([options.currentFolder '/output/' options.datasetName '/vb.mat'], 'vocabulary', 'allModes', 'distanceMatrices', 'modeProbs', 'options', '-append', '-v7.3');
   save([options.currentFolder '/output/' options.datasetName '/distributions.mat'], 'vocabularyDistributions', 'options');
   save([options.currentFolder '/output/' options.datasetName '/mainGraph.mat'], 'mainGraph', '-v7.3'); 
   
   %% Visualize images and vocabulary.
   if ~isempty(vocabLevel) && options.debug
     display('........ Visualizing previous levels...');
     [allNodeInstances, representativeNodes] =visualizeLevel( vocabLevel, vocabulary, graphLevel, firstLevelActivations, leafNodes, leafNodeCoords, levelItr, numel(vocabulary{1}), numel(vocabulary{levelItr-1}), options);
     display('........ Visualizing realizations on images...');
     if ~isempty(vocabLevel)
          matlabpool close;
          % Visualize realizations on images, and crop relevant
          % image parts.
          visualizeImages( fileList, vocabLevel, graphLevel, representativeNodes, allNodeInstances, leafNodes, leafNodeCoords, levelItr, options, 'train' );
          visualizeCroppedImgs( vocabLevel, representativeNodes, levelItr, options);

          % Backproject parts to the images.
          display('........ Imagining parts and their instances! This can take a while...');
          if options.vis.printTrainRealizations && levelItr >= 3
              try
                 projectTrainingImages(fileList, vocabularyDistributions, mainGraph, levelItr, options);
              catch
                 display('Visualization error'); 
              end
          end
          matlabpool('open', options.numberOfThreads);
     end
     printCloseFilters(eucDistanceMatrix, representativeNodes, levelItr, options); 
     clear allNodeInstances representativeNodes;
     java.lang.System.gc();
   end
   
   %% If we've reached max number of layers, don't keep going forward.
   if levelItr == options.maxLevels || stopFlag
       save([options.currentFolder '/output/' options.datasetName '/mainGraph.mat'], 'mainGraph', '-v7.3'); 
       return;
   end

   % Open/close matlabpool to save memory.
   matlabpool close;
   matlabpool('open', options.numberOfThreads);

   %% Step 2.6: If no new edges found, kill program.
   newEdgesAvailable = ~isempty(cat(1, mainGraph{levelItr}.adjInfo));
   if ~newEdgesAvailable
       vocabulary = vocabulary(1:(levelItr),:);
       vocabularyDistributions = vocabularyDistributions(1:(levelItr),:);
       allModes = allModes(1:(levelItr), :);
       mainGraph = mainGraph(1:(levelItr),:);
       distanceMatrices = distanceMatrices(1:(levelItr),:);
       modeProbs = modeProbs(1:(levelItr),:);
   end
   
   % clearing all, and then returning back to computation.
   scriptName = [options.currentFolder '/' options.datasetName '.sh'];
   restartFlag = options.restartFlag;
   if restartFlag
        % Finally, we save the workspace and call the next iteration.
         disp(['Saving layer ' num2str(levelItr) ' workspace.']);
         save([pwd '/Workspace.mat'], '-v7.3');
         copyfile([pwd '/Workspace.mat'], [options.currentFolder '/output/' options.datasetName '/Workspace' num2str(levelItr) '.mat']);
         clearvars -except restartFlag scriptName datasetName
         system(scriptName);
         exit
   else
         learnVocabularyLevel();
   end