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
  if exist('Workspace.mat', 'file')
       disp('Loading workspace from the previous layer.');
       loadWorkspace();
  end
   
    % Open threads for parallel processing.
    if options.parallelProcessing && usejava('jvm')
        if exist('parpool', 'file')
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
    
   % Increase level itr.
   levelItr = levelItr + 1; %#ok<*NASGU>
   
   % If this is a fresh Matlab, modify paths.
   % ========== PATH FOLDER ADDITION ==========
   if isempty(strfind(path, [options.currentFolder '/demo']))
       w = warning('off', 'all');
       addpath(genpath([options.currentFolder '/demo']));
       addpath(genpath([options.currentFolder '/graphTools']));
       addpath(genpath([options.currentFolder '/vocabLearning']));
       addpath(genpath([options.currentFolder '/inference']));
       addpath(genpath([options.currentFolder '/categorization']));
       warning(w);
   end
   
   if options.restartFlag
        % Try to restore break points.
        try %#ok<TRYNC>
            dbstop(bInfo);
        end
   end
   
   %% Before we start, we save relevant info regarding last layer so we can use it forming next layer.
    previousLevelImageIds = cat(1, graphLevel.imageId);
    previousLevelLeafNodes = {graphLevel.leafNodes};
    previousLevelPrecisePositions = cat(1, graphLevel.precisePosition);
    previousLevelPositions = cat(1, graphLevel.position);
    prevRealLabelIds = [graphLevel.realLabelId]';
    prevActivations = [graphLevel.activation]';
    
    %% Changing our learning architecture! 
    tic;
    [vocabLevel, graphLevel, vocabLevelDistributions] = discoverJointSubs(graphLevel, level1Coords, categoryArrIdx, options, levelItr-1);
    toc;
    
   %% Step 2.1: Run knowledge discovery to learn frequent compositions.
%   [vocabLevel, graphLevel, isSupervisedSelectionRunning, previousAccuracy] = discoverSubs(vocabLevel, graphLevel, level1Coords,...
%       options, levelItr-1, supervisedSelectionFlag, isSupervisedSelectionRunning, previousAccuracy); %#ok<*ASGLU,*NODEF>
%   if usejava('jvm')
%     java.lang.System.gc();
%   end

   % Open/close matlabpool to save memory.
%   matlabpool close;
%   matlabpool('open', options.numberOfThreads);

   %% If no new subs have been found, finish processing.
   if isempty(vocabLevel)
      % Write previous level's appearances to the output folder.
      vocabulary = vocabulary(1:(levelItr-1),:);
      vocabularyDistributions = vocabularyDistributions(1:(levelItr-1),:);
      distanceMatrices = distanceMatrices(1:(levelItr-1),:);
      return; 
   end

   %% Assign realizations R of next graph level (l+1), and fill in their bookkeeping info.
   graphLevel = fillBasicInfo(previousLevelImageIds, previousLevelLeafNodes, previousLevelPrecisePositions, graphLevel);
   if usejava('jvm')
     java.lang.System.gc();
   end

   %% Calculate statistics from this graph.
   display('........ Before we apply inhibition, estimating statistics..');
   [avgShareability, avgCoverage] = saveStats(vocabLevel, graphLevel, leafNodeCoords, maxCoverageVals, numberOfImages, options, 'preInhibition', levelItr);
   display(['........ Average Coverage: ' num2str(avgCoverage) ', average shareability of compositions: ' num2str(avgShareability) ' percent.']); 

   %% If category level is reached, we reduce the number of desired nodes substantially. 
   display('........ Checking for image coverages for category level..');
   imageIds = [graphLevel.imageId];
   [options, nodeCoverages, imageCoverages] = markCategoryLevel(graphLevel, level1Coords, avgCoverage, levelItr, options);

   %% In order to do proper visualization, we learn precise positionings of children for every vocabulary node.
   display('........ Learning sub-part label and position distributions.');
   [vocabLevel, vocabLevelDistributions] = learnChildDistributions(vocabLevel, vocabLevelDistributions, graphLevel, prevRealLabelIds, single(previousLevelPositions), previousLevelPrecisePositions, levelItr, options);

   if usejava('jvm')
     java.lang.System.gc();
   end

   %% Calculate activations for every part realization.
   display('........ Calculating activations..');
   [vocabLevel, vocabLevelDistributions, graphLevel] = calculateActivations(vocabLevel, vocabLevelDistributions, graphLevel, prevActivations, single(previousLevelPositions), levelItr, options);
   if usejava('jvm')
     java.lang.System.gc();
   end
   
   %% Post-process graphLevel, vocabularyLevel to remove non-existent parts from vocabLevel.
   % In addition, we re-assign the node ids in graphLevel.
  display('........ Calculating distance matrix among the vocabulary nodes (in parallel)..');
  vocabularyDistributions{levelItr} = vocabLevelDistributions;   
  [vocabLevel, vocabularyDistributions, graphLevel, eucDistanceMatrix] = postProcessParts(vocabLevel, graphLevel, vocabularyDistributions, levelItr, options);    

  % Visualize learned OR Nodes.
  visualizeORNodes( vocabLevel, levelItr, options);
  
  % Update distance matrices.
  newDistanceMatrix = inf(size(eucDistanceMatrix,1), size(eucDistanceMatrix,1), 'single');
  newDistanceMatrix(1:(size(eucDistanceMatrix,1)+1):size(eucDistanceMatrix,1)*size(eucDistanceMatrix,1)) = 0;
  distanceMatrices{levelItr} = newDistanceMatrix;
  if usejava('jvm')
     java.lang.System.gc();
  end

   %% We process graphLevel's labelIds to reflect updated labels (OR Node Labels).
   vocabNodeLabels = [vocabLevel.label];
   updatedLabelIds = num2cell(vocabNodeLabels([graphLevel.labelId]));
   [graphLevel.labelId] = deal(updatedLabelIds{:});
   clear updatedLabelIds;
   
   %% As an initial stage of inhibition, we downsample the responses, 
   % and apply max pooling. Please note that turning this off also means no
   % reduction in resolution for higher layers, therefore non-growing
   % receptive fields.
   if ~ismember(levelItr, options.noPoolingLayers)
        positions = cat(1, graphLevel.position);
        positions = poolPositions(positions, options.poolDim);
        newPositions = cell(numel(graphLevel),1);
        for itr = 1:numel(graphLevel)
             newPositions{itr} = positions(itr, :);
        end
        [graphLevel.position] = deal(newPositions{:});
        level1Coords(:, 4:5) = poolPositions(level1Coords(:, 4:5), options.poolDim);
   end
   display(['........ Applying pooling on ' num2str(numel(graphLevel)) ' realizations belonging to ' num2str(max([vocabLevel.label])) ' compositions.']);
   graphLevel = applyPooling(graphLevel, options.poolFlag, options.labeledPooling);
   display(['........ After pooling, we have ' num2str(numel(graphLevel)) ' realizations of ' num2str(numel(unique([graphLevel.labelId]))) ' compositions.']);
   if usejava('jvm')
     java.lang.System.gc();
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
       display(['........ Inhibition finished. We have ' num2str(numel(graphLevel)) ' realizations of ' num2str(numel(unique([graphLevel.labelId]))) ' compositions.']);
   else
       display('........ No inhibition applied..');
   end
   
   %% Post-inhibition check: If all realizations of a vocabulary node are missing, we can move on.
   % Assign new labels of the remaining realizations.
   [graphLevel, vocabLevel, vocabularyDistributions] = updateDataStructures(graphLevel, vocabLevel, vocabularyDistributions, levelItr, options);
   
   %% Calculate statistics from this graph.
   display('........ Estimating post-inhibition statistics..');
   [avgShareability, avgCoverage] = saveStats(vocabLevel, graphLevel, leafNodeCoords, maxCoverageVals, numberOfImages, options, 'postInhibition', levelItr);

   % display debugging info.
   display(['........ Remaining: ' num2str(numel(graphLevel)) ' realizations belonging to ' num2str(max([vocabLevel.label])) ' compositions.']);
   display(['........ Average Coverage: ' num2str(avgCoverage) ', average shareability of compositions: ' num2str(avgShareability) ' percent.']); 

   %% Step 2.5: Create the parent relationships between current level and previous level.
   vocabulary(levelItr) = {vocabLevel};

   if usejava('jvm')
     java.lang.System.gc();
   end
   
   %% Print distance matrices.
   if ~isempty(newDistanceMatrix)
      imwrite(eucDistanceMatrix, [options.currentFolder '/debug/' options.datasetName '/level' num2str(levelItr) '_dist.png']);
      imwrite(newDistanceMatrix, [options.currentFolder '/debug/' options.datasetName '/level' num2str(levelItr) '_orNodes.png']);
   end
   
   %% At this step, we perform imagination of the parts at this layer (for later use).
   vocabularyDistributions = projectVocabulary( vocabularyDistributions, options);

   %% If we've reached max number of layers, don't keep going forward.
   if levelItr == options.maxLevels || options.stopFlag
       vocabulary = vocabulary(1:(levelItr),:);
       vocabularyDistributions = vocabularyDistributions(1:(levelItr),:);
       distanceMatrices = distanceMatrices(1:(levelItr),:);
   end
   
   %% We're exporting output here. This helps us to perform 
   % Export realizations into easily-readable arrays.
   display('Writing output to files.');
   [exportArr, activationArr, pooledPositions, allLeafNodes] = exportRealizations(graphLevel, levelItr);
   exportArr = cat(1, preExportArr, exportArr);
   preExportArr = exportArr;
   allLeafNodes = cat(1, preLeafNodes, allLeafNodes);
   preLeafNodes = allLeafNodes;
   pooledPositions = cat(1, prePooledPositions, pooledPositions);
   prePooledPositions = pooledPositions;
   activationArr = cat(1, preActivationArr, activationArr);
   preActivationArr = activationArr;
   save([options.currentFolder '/output/' options.datasetName '/export.mat'], 'exportArr', 'pooledPositions', 'activationArr', '-append'); 
   save([options.currentFolder '/output/' options.datasetName '/leafNodes.mat'], 'leafNodes');
   clear activationArr precisePositions;

   % Print everything to files.
   save([options.currentFolder '/output/' options.datasetName '/vb.mat'], 'vocabulary', 'distanceMatrices', 'options', '-append');
   save([options.currentFolder '/output/' options.datasetName '/distributions.mat'], 'vocabularyDistributions', 'options', '-v7');
   
   %% Visualize images and vocabulary.
   if ~isempty(vocabLevel) && options.debug
     display('........ Visualizing previous levels...');
     [allNodeInstances, representativeNodes] = visualizeLevel( vocabLevel, graphLevel, firstLevelActivations, leafNodes, leafNodeCoords, levelItr, numel(vocabulary{1}), numel(vocabulary{levelItr-1}), options);
     display('........ Visualizing realizations on images...');
     if ~isempty(vocabLevel)
 %         matlabpool close;
          % Visualize realizations on images, and crop relevant
          % image parts.
          visualizeImages( fileList, graphLevel, representativeNodes, allNodeInstances, leafNodes, leafNodeCoords, levelItr, options, 'train' );
          visualizeCroppedImgs( vocabLevel, representativeNodes, levelItr, options);

          % Backproject parts to the images.
          display('........ Imagining parts and their instances! This can take a while...');
          if options.vis.printTrainRealizations
              try
                 projectTrainingImages(fileList, vocabularyDistributions, exportArr, levelItr, options);
              catch
                 display('Visualization error'); 
              end
          end
%          matlabpool('open', options.numberOfThreads);
     end
%     printCloseFilters(eucDistanceMatrix, representativeNodes, levelItr, options); 
     clear allNodeInstances representativeNodes;
     if usejava('jvm')
        java.lang.System.gc();
     end
   end
   
   % Open/close matlabpool to save memory.
%   matlabpool close;
%   matlabpool('open', options.numberOfThreads);
   
   % clearing all, and then returning back to computation.
   scriptName = [options.currentFolder '/' options.datasetName '.sh'];
   restartFlag = options.restartFlag;
   
   %% If we've reached max number of layers, don't keep going forward.
   if levelItr == options.maxLevels || options.stopFlag
       return;
   elseif numel(fileList) >= options.saveWorkspaceImageCount
        % Save workspace for later operation.
        disp(['Saving layer ' num2str(levelItr) ' workspace.']);
        saveWorkspace();
   end
   
   % Restart program.
   if restartFlag
        bInfo = dbstatus;
        % Finally, we save the workspace and call the next iteration.
        clearvars -except restartFlag scriptName datasetName
        system(scriptName);
        exit
   else
         learnVocabularyLevel();
   end