%> Name: discoverJointSubs
%>
%> Description: This function is created in order to discover n-groups of
%> substructures in high dimensional space. Instead of creating duplets and
%> then mining groups of duplets, this method tries to find groups that
%> satisfy a multiple-criterion objective function, and merges duplet
%> learning with substructure discovery. 
%>
%> @param vocabLevel Parts.
%> @param graphLevel Part realizations.
%> @param level1Coords Layer 1 nodes' coordinates.
%> @param options Program options.
%> @param levelItr Current level id.
%> 
%> @retval vocabLevel Learned parts.
%> @retval graphLevel Learned realizations.
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 24.11.2016
function [newVocabLevel, newGraphLevel, newVocabLevelDistributions] = discoverJointSubs(graphLevel, level1Coords, categoryArrIdx, options, levelItr)
     %% Create input representation.
     display('Discovering new parts in parallel..');
     prevLevel = cat(2, cat(1, graphLevel.labelId), cat(1, graphLevel.position), cat(1, graphLevel.imageId));
     prevLevelPositions = single(cat(1, graphLevel.position));
     leafNodes = {graphLevel.leafNodes};
     rfSize = options.receptiveFieldSizes(levelItr);
     nodeCoverageThr = options.missingNodeThr;
     numberOfSelectedSubs = options.reconstruction.numberOfReconstructiveSubs(levelItr);
     selectionOptions.numberOfSelectedSubs = numberOfSelectedSubs;
     selectionOptions.unsupervisedWeight = 1;
     selectionOptions.supervisedWeight = 0;
     selectionOptions.reconstructivePartSelection = options.reconstructivePartSelection;
     selectionOptions.discriminativePartSelection = options.discriminativePartSelection;
     selectionOptions.leafNodeCover = options.leafNodeCover;
     minSize = 2;
     maxClusters = 10;
     
     newFolder = [options.debugFolder '/level' num2str(levelItr+1) '/candidateParts'];
     if ~exist(newFolder, 'dir')
          mkdir(newFolder);
     end
     w = warning('off', 'all');
     
     %% Obtain all cliques from the data (may be missing some parts if there are too many combinations)
     display(['Obtaining cliques by beam search from ' num2str(numel(unique(prevLevel(:,1)))) ' parts and ' num2str(size(prevLevel,1)) ' instances...']);
     [allCliques, labeledCliques] = getCombinations(prevLevel, rfSize, minSize, options.subdue.maxSize, options.subdue.beam, options.subdue.beamReductionRate, options.subdue.nsubs);

     % Obtain candidate parts.
     [candidateParts, firstInstances, idx] = unique(labeledCliques, 'rows', 'R2012a');
     clear labeledCliques;
     firstInstances = cat(1, firstInstances, numel(idx)+1);
     
     %% Calculate instance positions!
     % This process removes parts that have consistently low coverage!
     prevLevelX = prevLevelPositions(:,1);
     prevLevelY = prevLevelPositions(:,2);
     
    allInstancePositions = cell(size(candidateParts,1),1);
    parfor partItr = 1:size(candidateParts, 1)             
          % Get relevant instances.
          validSubPartCount = nnz(candidateParts(partItr,:));
          relevantCliques = allCliques(firstInstances(partItr):(firstInstances(partItr+1)-1),:);
          relevantCliques = relevantCliques(:, (end + (1 - validSubPartCount)):end);
          
          % Calculate centers
          numberOfCliques = size(relevantCliques,1);
          instancePositions = zeros(numberOfCliques,2, 'single');
          xvals = reshape(prevLevelX(relevantCliques), size(relevantCliques));
          instancePositions(:,1) = (min(xvals, [], 2) + max(xvals, [], 2)) / 2;
          yvals = reshape(prevLevelY(relevantCliques), size(relevantCliques));
          instancePositions(:,2) = (min(yvals, [], 2) + max(yvals, [], 2)) / 2;
          allInstancePositions{partItr} = round(instancePositions);
    end
     allInstancePositions = cat(1, allInstancePositions{:});
   
     %% Joing Statistical Learning to obtain representative parts.
     % Part selection and optimization is performed in this crucial step.
     % Instead of duplets, parts are learned jointly. 
     newCliques = cell(size(candidateParts,1),1);
     newLabels = cell(size(candidateParts,1),1);
     newInstancePositions = cell(size(candidateParts,1),1);
     newPartDistributions = cell(size(candidateParts,1),1);
     newPartThresholds = cell(size(candidateParts, 1),1);
     newChildrenArr = cell(size(candidateParts,1),1);
 
     display(['Discovering new parts by learning multi-modal distributions... We have ' num2str(size(candidateParts,1)) ' different types of cliques.']);
     parfor partItr = 1:size(candidateParts,1)
%     for partItr = 1:200
          w2 = warning('off', 'all');
          instancePositions = allInstancePositions(firstInstances(partItr):(firstInstances(partItr+1)-1),:);
          %% Create cliques and their joint position space.
          [jointPositions, savedCliques, uniqueSubParts, numberOfSubParts] = getCliques(candidateParts, allCliques, firstInstances, prevLevelPositions, allInstancePositions, rfSize, partItr);
          %% Calculate joint statistics. 
          [refinedClusterSamples, refinedCliques, refinedClusters, refinedInstancePositions, shownSamples, shownClusters, clusterDistributions, clusterThresholds, numberOfClusters, partClusters] = learnJointStats(jointPositions, instancePositions, savedCliques, rfSize, maxClusters);
          newCliques{partItr} = refinedCliques;
          newLabels{partItr} = int32(refinedClusters);
          newInstancePositions{partItr} = refinedInstancePositions;
          newPartDistributions{partItr} = clusterDistributions;
          newPartThresholds{partItr} = clusterThresholds;
          
          %% Visualize parts.
%          visualizeJointParts(shownSamples, shownClusters, refinedClusterSamples, refinedClusters, jointPositions, partClusters, numberOfSubParts, partItr, rfSize, numberOfClusters, uniqueSubParts, newFolder);
          children = candidateParts(partItr,:);
          newChildrenArr{partItr} = repmat(children, numel(unique(refinedClusters)),1);
          warning(w2);
     end
     
     % Update part labels.
     offset = 0;
     for partItr = 1:size(candidateParts,1)
          newLabels{partItr} = newLabels{partItr} + offset;
          tempArr = newLabels{partItr};
          if ~isempty(tempArr)
               offset = offset + numel(fastsortedunique(sort(newLabels{partItr})));
          end
     end
     
     % Save all info and proceed to selection.
     newCliques = cat(1, newCliques{:});
     newLabels = cat(1, newLabels{:});
     newInstancePositions = cat(1, newInstancePositions{:});
     newPartDistributions = cat(1, newPartDistributions{:});
     newPartThresholds = cat(1, newPartThresholds{:});
     newChildrenArr = cat(1, newChildrenArr{:});
     newImageIds = prevLevel(newCliques(:, end),4);
     newCategoryLabels = categoryArrIdx(newImageIds);
     newMetaData = cat(2, newLabels, newImageIds, newCategoryLabels);
     clear newLabels newImageIds newCategoryLabels;
     
     %% Before the selection process, we eliminate parts which have low coverage. 
     % They are removed from the part selection process.
     if nodeCoverageThr > 0
          display('Removing low coverage parts...');
          [newCliques, newInstancePositions, newMetaData, newChildrenArr, newPartDistributions, newPartThresholds] = removeLowCoverageParts(newCliques, newInstancePositions, newMetaData, newChildrenArr, newPartDistributions, newPartThresholds, prevLevel, leafNodes, level1Coords, nodeCoverageThr, rfSize);
     end
     
     %% Candidate parts are discovered. Now, we perform part selection process.
     str = 'Part selection is being performed. Part selection modes:';
     if selectionOptions.reconstructivePartSelection
          str = [str ' [Reconstructive]'];
     end
     if selectionOptions.discriminativePartSelection
          str = [str ' [Discriminative]'];
     end
     disp(str);
     selectedParts = selectJointParts(newCliques, newMetaData, leafNodes, level1Coords, categoryArrIdx, selectionOptions); % newCliques, sub positions, level 1 positions, category labels, image ids, Remaining level 1 nodes
     
     %% Create new data structures from selected parts. 
     idx = ismembc(newMetaData(:,1), int32(selectedParts));
     display(['Discovered ' num2str(numel(selectedParts)) ' parts with ' num2str(nnz(idx)) ' instances. Writing output and exit.']);
     newCliques = newCliques(idx, :);
     newInstancePositions = newInstancePositions(idx,:);
     newMetaData = newMetaData(idx, :);
     [~, ~, updatedLabels] = unique(newMetaData(:,1));
     numberOfParts = numel(selectedParts);
     numberOfInstances = nnz(idx);
     assignedLabels = num2cell(int32(updatedLabels));
     assignedVocabLabels = num2cell(int32((1:numberOfParts)'));
     childrenLabels = newChildrenArr(selectedParts,:);
     newPartDistributions = newPartDistributions(selectedParts, :);
     newPartThresholds = newPartThresholds(selectedParts, :);
     newPartThresholds = num2cell(newPartThresholds);
     childrenLabels = mat2cell(childrenLabels, ones(numberOfParts,1), size(childrenLabels,2));
     childrenLabels = cellfun(@(x) x(x>0), childrenLabels, 'UniformOutput', false);
     
     % Allocate space for new stuff.
     newVocabLevel(numberOfParts) = VocabNode();
     newVocabLevelDistributions(numberOfParts) = NodeDistribution();
     newGraphLevel(numberOfInstances) = GraphNode();
     
     % Assign relevant values.
     childrenArr = cell(numel(newGraphLevel),1);
     positionArr = cell(numel(newGraphLevel),1);
     for itr = 1:numel(newGraphLevel)
          tmp = newCliques(itr, :);
          tmp = tmp(tmp>0);
          childrenArr{itr} = tmp;
          positionArr{itr} = int32(newInstancePositions(itr, :));
     end
     [newGraphLevel.children] = deal(childrenArr{:});
     [newGraphLevel.labelId] = deal(assignedLabels{:});
     [newGraphLevel.realLabelId] = deal(assignedLabels{:});
     [newGraphLevel.position] = deal(positionArr{:});
     [newGraphLevel.sign] = deal(1);
     [newGraphLevel.activation] = deal(single(1));
     
     [newVocabLevel.label] = deal(assignedVocabLabels{:});
     [newVocabLevel.children] = deal(childrenLabels{:});
     
     [newVocabLevelDistributions.childrenPosDistributions] = deal(newPartDistributions{:});
     [newVocabLevelDistributions.distanceThr] = deal(newPartThresholds{:});
     
     warning(w);
end