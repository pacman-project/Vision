%> Name: runSubdue
%>
%> Description: Inference on graphLevel with given subs in vocabLevel. Each
%> composition is searched for in graphLevel.
%>
%> @param vocabLevel Input vocabulary level. Compositions in this vocabulary 
%> level are detected in graphLevel.
%> @param graphLevel The current object graphs' level. The graphLevel's
%> nodes are sorted first by their imageId, then labelId.
%> @param options Program options.
%>
%> @retval graphLevel The graph level consisting of discovered nodes.
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 05.02.2014
function [exportArr, activationArr] = inferSubs(fileName, img, vocabulary, vocabularyDistributions, allModes, modeProbs, nodes, nodeActivations, options)
    % Read data into helper data structures.
    RFSize = options.receptiveFieldSize;
    smallHalfMatrixSize = (options.smallReceptiveFieldSize+1)/2;
    missingPartLog = log(0.00001);
    maxShareability = options.maxShareability;
    missingPartAllowed = false;
    poolFlag = true;
    maxPosProb = 0.1;
    halfMatrixSize = (options.receptiveFieldSize+1)/2;
    load([pwd '/filters/optimizationFilters.mat'], 'visFilters');
    outputImgFolder  = [options.currentFolder '/output/' options.datasetName '/reconstruction/test/' fileName];
    filterIds = round(((180/numel(options.filters)) * (0:(numel(options.filters)-1))) / (180/size(visFilters,3)))' + 1;
    
    exportArr = [];
    if isempty(nodes) || isempty(vocabulary)
        return;
    end
    
    % Output data structures.
    allNodes = cell(numel(vocabulary),1);
    allActivations = cell(numel(vocabulary),1);
    allPrecisePositions = cell(numel(vocabulary),1);
    
    % Fill layer 1 data structures.
    precisePositions = single(nodes(:,4:5));
    leafPrecisePositions = precisePositions;
    leafNodes = num2cell((1:size(nodes,1))');
    
    allNodes(1) = {nodes(:,1:3)};
    allActivations(1) = {zeros(size(nodes,1),1,'single')};
    allLeafNodes{1} = leafNodes;
    allPrecisePositions{1} = precisePositions;
    prevVocabORLabels = cat(1, vocabulary{1}.label);
    
    if strcmp(options.filterType, 'gabor')
         stride = options.gabor.stride;
    else
         stride = options.auto.stride;
    end
    poolDim = options.poolDim;
    halfRFSize = ceil(RFSize/2);
    w = warning('off', 'all');
    
%     rfRadius = halfRFSize;
%     minRFRadius = 3;
%     %% Create a network of first layer nodes, and their links for upper layers.
%     firstLayerDistances = squareform(pdist(single(nodes(:,2:3))));
%     allAdjacentNodes = cell(size(nodes,1), 1);
%     for nodeItr = 1:size(nodes,1)
%          distances = firstLayerDistances(nodeItr, :);
%          adjacentNodes = distances < halfRFSize & distances >= minRFRadius;
%          adjacentNodes(nodeItr) = 0;
%          adjacentNodes = find(adjacentNodes);
%     end
     
    %% Start discovery process.
    for vocabLevelItr = 2:min(numel(vocabulary), 10)    
         rfSize = getRFSize(options, vocabLevelItr);
         halfRealRFSize = round(rfSize(1)/2);
         realRFSize = halfRealRFSize * 2 + 1;
         modes = allModes{vocabLevelItr-1};
         curModeProbs = modeProbs{vocabLevelItr-1};
         [rankModes, ~, ~] = unique(modes(:,1:2), 'rows', 'R2012a');
         rankModeIds = sparse(double(rankModes(:,1)), double(rankModes(:,2)), 1:size(rankModes,1));
         
         if ismember(vocabLevelItr-1, options.smallRFLayers)
             rfRadius = smallHalfMatrixSize;
         else
             rfRadius = halfRFSize;
         end
         
         % Obtain edge novelty threshold.
         edgeNoveltyThr = max(options.minEdgeNoveltyThr, options.edgeNoveltyThr - options.edgeNoveltyThrRate * max(0, (vocabLevelItr-3)));
         edgeShareabilityThr = 1 - edgeNoveltyThr;
         allowedSharedLeafNodes = cellfun(@(x) numel(x), leafNodes) * edgeShareabilityThr;
         
        % If no pooling has been performed at this layer, and previous layer
        % has small RF, we have a minimum RF limit. 
        if ismember(vocabLevelItr-1, options.noPoolingLayers) && ismember(vocabLevelItr - 2, options.smallRFLayers) && vocabLevelItr < 8
            minRFRadius = floor(halfMatrixSize/2);
        else
%            minRFRadius = 1;
            minRFRadius = max(1, min(3, 5-vocabLevelItr));
        end
         
        % Calculate pool factor.
        poolFactor = nnz(~ismembc(2:vocabLevelItr, options.noPoolingLayers));
         
        %% Match subs from vocabLevel to their instance in graphLevel.
        % Obtain relevant data structures for the level.
        vocabLevel = vocabulary{vocabLevelItr};
        vocabLevelDistributions = vocabularyDistributions{vocabLevelItr};
        if isempty(vocabLevel)
             break;
        end
        
        % Allocate space for this layer's realizations.
        layerNodes = cell(numel(vocabLevel),1);
        layerActivations = layerNodes;
        layerLeafNodes  = layerNodes;
        layerPrecisePositions = layerNodes;
        
        % Obtain OR node labels.
        nodeORLabels = prevVocabORLabels(nodes(:,1));
        
        % Calculate pairwise distances.
        if size(nodes,1) == 1
             distances = Inf;
        else
             distances = squareform(pdist(single(nodes(:,2:3))));
             distances(1:(size(distances,1)+1):size(distances,1)^2) = Inf;
        end
        
        %% Put statistics in cell arrays for faster reference.
         posArr = {vocabLevelDistributions.childrenPosDistributionProbs};
         minPosLogProbs = cat(1, vocabLevelDistributions.minPosActivationLog);
         minLogProbs = cat(1, vocabLevel.minActivationLog);
         
        %% Starting inference from layer 2. 
         for vocabItr = 1:numel(vocabLevel)
               % Match the center first.
               vocabNode = vocabLevel(vocabItr);
               vocabNodeChildren = vocabNode.children;
               centerId = vocabNodeChildren(1);
               validInstances = nodeORLabels == centerId;

%                if numel(vocabNodeChildren) == 1
%                   continue;
%                end

               % Obtain edge matrices.
               if numel(vocabNodeChildren) > 1
                   secChildren = vocabNodeChildren(2:end)';
                   matrixIdx = zeros(numel(secChildren),1);
                   for itr = 1:numel(matrixIdx)
                        matrixIdx(itr) = rankModeIds(vocabNodeChildren(1), secChildren(itr));
                   end
                   edgeMatrices = curModeProbs(:,:,matrixIdx);
               end

               % If no centers have matched, continue.
               if nnz(validInstances) == 0
                 continue;
               end

               % Allocate space to hold the instances and their activations.
               centerNodes = find(validInstances);
               numberOfCenterNodes = numel(centerNodes);
               instanceChildren = cell(numberOfCenterNodes,numel(vocabNodeChildren));

               %% Iteratively match instances and get all node combinations.
               validCombinations = ones(numberOfCenterNodes,1) > 0;
               for centerItr = 1:numberOfCenterNodes
                  curCenterNode = centerNodes(centerItr);
                  instanceChildren(centerItr, 1) = {curCenterNode};
                  validCols = ones(numel(vocabNodeChildren),1) > 0;
                  
                  % Search for secondary nodes.
                  for childItr = 2:numel(vocabNodeChildren)
                        tempNodes = find(nodeORLabels == vocabNodeChildren(childItr) & distances(:, curCenterNode) < rfRadius & distances(:, curCenterNode) >= minRFRadius);
                        
                        % Remove repetitions.
                        tempNodes = tempNodes(tempNodes ~= curCenterNode);

                        if ~isempty(tempNodes)
                             % Fine-tune the system by filtering temp nodes
                             % more.
                             curEdgeType = vocabNode.adjInfo(childItr-1, 3);
                             relativeCoords = nodes(tempNodes,2:3) + halfRFSize;
                             relativeCoords(:,1) = relativeCoords(:,1) - nodes(curCenterNode, 2);
                             relativeCoords(:,2) = relativeCoords(:,2) - nodes(curCenterNode, 3);
                             linIdx = relativeCoords(:,1) + (relativeCoords(:,2)-1)*RFSize;
                             curEdgeMatrix = edgeMatrices(:,:,childItr-1);

                             % Save edge types.
                             edgeTypes = curEdgeMatrix(linIdx);
                             tempNodes = tempNodes(edgeTypes == curEdgeType);
                        end
                        
                        % Edge novelty threshold applied here.
                        if ~isempty(tempNodes) && vocabLevelItr > 2
                             if edgeShareabilityThr < 1
                                  commonLeafCounts = cellfun(@(x) sum(ismembc(x, leafNodes{curCenterNode})), leafNodes(tempNodes));
                                  novelNodes = (commonLeafCounts <= allowedSharedLeafNodes(tempNodes)) & ...
                                      (commonLeafCounts <= allowedSharedLeafNodes(curCenterNode));
                                  tempNodes = tempNodes(novelNodes);
                             end
                        end

                        % If tempNodes is not empty, save it.
                        if isempty(tempNodes)
                 %            if missingPartAllowed
                  %                tempNodes = Inf;
                  %                validCols(childItr) = 0;
                  %           else
                                  validCombinations(centerItr) = 0;
                                  break;
                  %           end
                        end
                        instanceChildren(centerItr, childItr) = {tempNodes};
                  end
               end
               
               % Create node combinations. We have different solutions
               % based on the number of repetitions.
               instanceChildren = instanceChildren(validCombinations, :);
               cellCounts = cellfun(@(x) numel(x), instanceChildren);
               instanceCombinations = cell(size(instanceChildren,1),1);
               for instanceItr = 1:size(instanceChildren,1)
                    multipleVals = find(cellCounts(instanceItr,:) > 1);
                    numberOfMultipleVals = numel(multipleVals);
                    % No repetitions? Add them side by side.
                    if numberOfMultipleVals == 0
                         curInstanceCombinations = cat(2, instanceChildren{instanceItr, :});
                    elseif numberOfMultipleVals == 1
                         % One repetition is a simple case where one column
                         % has multiple choices, and others are repeated.
                         curInstanceCombinations = zeros(cellCounts(instanceItr, multipleVals), size(instanceChildren,2));
                         curInstanceCombinations(:, multipleVals) = instanceChildren{instanceItr, multipleVals};
                         otherCols = 1:size(instanceChildren,2);
                         otherCols = otherCols(otherCols ~= multipleVals);
                         for itr = otherCols
                              curInstanceCombinations(:,itr) = instanceChildren{instanceItr, itr};
                         end
                    else
                         % Tricky case with multiple columns with multiple
                         % choices. Use allcomb.
                         curInstanceCombinations = allcomb(instanceChildren{instanceItr, :});
                    end
                    instanceCombinations{instanceItr} = curInstanceCombinations;
               end
               
               % Find all node combinations (not containing node
               % repetitions). Stale, will be used when missing nodes are
               % introduced.
%                curInstanceCombinations = allcomb(instanceChildren{centerItr,validCols});
%                sortedmatrix = sort(curInstanceCombinations,2);
%                validInstanceIdx = all(diff(sortedmatrix,[],2)~=0,2);
%                if nnz(validCols) ~= numel(validCols)
%                   curInstanceCombinations = allcomb(instanceChildren{centerItr,:});
%                end
%                curInstanceCombinations = curInstanceCombinations(validInstanceIdx,:);
%                instanceCombinations{centerItr} = curInstanceCombinations;
               
               % We have the instances for this sub. Save the info.
               instanceCombinations = cat(1, instanceCombinations{:});
               numberOfCombinations = size(instanceCombinations,1);
               numberOfChildren = size(instanceCombinations,2);
               
               if numberOfCombinations == 0
                    continue;
               end
               
              %% Check leaf node shareability, and save leaf nodes for every instance.
              instanceLeafNodes = cell(numberOfCombinations, 1);
              shareabilityArr = zeros(numberOfCombinations,1);
              for instanceItr = 1:numberOfCombinations
                   tempCombinations = instanceCombinations(instanceItr,:);
                   
                   tempArr = sort(cat(1, leafNodes{tempCombinations(~isinf(tempCombinations))}));
                   tempLeafNodes = fastsortedunique(tempArr);
%                    
%                    if maxShareability < 1
%                         numberOfUniqueLeafNodes = numel(tempLeafNodes);
%                         numberOfRepetitions = numel(tempArr) - numel(tempLeafNodes);
% 
%                         % Determine if this instance is any good.
%                         repetitionRatio = numberOfRepetitions / numberOfUniqueLeafNodes;
%                         if repetitionRatio > maxShareability * (numberOfChildren-1)
%                              shareabilityArr(instanceItr) = 1;
%                         elseif repetitionRatio > maxShareability
%                              tempArr = [false; diff(tempArr)== 0];
%                              tempArr = find(tempArr);
%                              repeatedNodes = numel(tempArr);
%                              for itr = 2:numel(tempArr)
%                                    if tempArr(itr) == tempArr(itr-1)+1
%                                         repeatedNodes = repeatedNodes - 1;
%                                    end
%                              end
%                              shareabilityArr(instanceItr) = repeatedNodes / numberOfUniqueLeafNodes;
%                          end
%                    end
                   instanceLeafNodes{instanceItr} = tempLeafNodes;
              end
              
              % Filter invalid instances by shareability threshold.
              validInstances = shareabilityArr <= maxShareability;
              if nnz(validInstances) == 0
                   continue;
              end
              instanceCombinations = instanceCombinations(validInstances,:);
              instanceLeafNodes = instanceLeafNodes(validInstances,:);
              numberOfCombinations = size(instanceCombinations,1);
              numberOfChildren = size(instanceCombinations,2);
               
              %% Combinations are discovered. Now, we calculate an activation value for every combination.
              % First, we have to pick relevant distributions for each node. We
              % need combined relative positions first.
              % Calculate instance positions.
              instancePrecisePositions = zeros(numberOfCombinations,2, 'single');
              for instanceItr = 1:numberOfCombinations
                  positions = leafPrecisePositions(instanceLeafNodes{instanceItr},:);
                  instancePrecisePositions(instanceItr,:) = round(sum(positions,1)/size(positions,1));
              end
             
             %% Calculate activations!
             childrenProbs = posArr{vocabItr};
             
             childrenPosActivations = zeros(size(instanceCombinations),1, 'single');
             for childItr = 1:size(instanceCombinations,2)
                  relevantChildren = instanceCombinations(:, childItr);
                  relevantProbs = childrenProbs{childItr};

                  % Calculate position probabilities.
                  samples = (precisePositions(relevantChildren, :) - instancePrecisePositions) + halfRealRFSize + 1;
                  validIdx = samples > 0 & samples < realRFSize + 1;
                  validIdx = validIdx(:,1) & validIdx(:,2);
                  pointProbs = zeros(size(samples,1),1, 'single');
                  samplesIdx = samples(validIdx,1) + (samples(validIdx,2)-1)*realRFSize;
                  pointProbs(validIdx) = full(relevantProbs(samplesIdx));
                  childrenPosActivations(:, childItr) = log(pointProbs);
             end
             instancePosActivations = sum(childrenPosActivations,2) / numberOfChildren;
             
             % Take previous activations into account as well.
             prevChildrenActivations = nodeActivations(instanceCombinations);
             if ~isequal(size(prevChildrenActivations), size(instanceCombinations))
                 prevChildrenActivations = prevChildrenActivations';
             end
             if vocabLevelItr > 2
                   instanceActivations = sum([prevChildrenActivations, childrenPosActivations], 2) / (numberOfChildren * 2);
             else
                   instanceActivations = sum(childrenPosActivations,2) / numberOfChildren;
             end
             
             % Find valid activations by filtering.
             validCombinations = instancePosActivations >= minPosLogProbs(vocabItr) & instanceActivations >= minLogProbs(vocabItr);
             
              %% Filter the data structures by eliminating invalid (low probability) combinations.
             instanceActivations = instanceActivations(validCombinations,:);
             instancePrecisePositions = instancePrecisePositions(validCombinations,:);
             instanceLeafNodes = instanceLeafNodes(validCombinations,:);
             
             %% Create new nodes.
             newNodes = [ones(nnz(validCombinations), 1, 'int32') * vocabItr, int32(instancePrecisePositions)];

             %% Finally, write everything to the relevant arrays.
             layerNodes{vocabItr} = newNodes;
             layerActivations{vocabItr} = instanceActivations;
             layerLeafNodes{vocabItr}  = instanceLeafNodes;
             layerPrecisePositions{vocabItr} = instancePrecisePositions;
         end
         
         %% Processing of a layer is finished. We have to move on to the next layer here.
         % First, pooling.
         nodes = cat(1, layerNodes{:});
         nodeActivations = cat(1, layerActivations{:});
         leafNodes = cat(1, layerLeafNodes{:});
         precisePositions = cat(1, layerPrecisePositions{:});
         prevVocabORLabels = cat(1, vocabLevel.label);
         
         if isempty(nodes)
              break;
         end
         
          %% Perform likelihood-based inhibition
          pooledPositions = calculatePooledPositions(precisePositions, poolFactor, poolDim, stride );
         
          if isempty(nodes)
              break;
          end
         
          % Get unique nodes for each label, coords, activation triplet.
          combinedArr = double([-nodeActivations, double(prevVocabORLabels(nodes(:,1))), double(pooledPositions)]);

          % First, we order the nodes.
          [combinedArr, sortIdx] = sortrows(combinedArr);
          
          %% Perform max pooling and save data structures.
          if poolFlag
               [~, IA, ~] = unique(combinedArr(:,2:4), 'rows', 'stable');
          else
               IA = (1:size(combinedArr,1))';
          end

          % Save remaining data and move on to the next layer.
         validIdx = sort(sortIdx(IA));
         nodes = nodes(validIdx,:);
         nodeActivations = nodeActivations(validIdx,:);
         leafNodes = leafNodes(validIdx, :);
         precisePositions = precisePositions(validIdx, :);
         nodes(:,2:3) = calculatePooledPositions(precisePositions, poolFactor, poolDim, stride );
         
         %% For debugging purposes, we visualize the nodes.
         if options.testDebug
              experts = projectNode([nodes(:,1), precisePositions, ones(size(nodes,1),1)*vocabLevelItr], vocabularyDistributions, 'modal', options);
              experts(:,1) = filterIds(experts(:,1));
              modalImg = obtainPoE(experts, [], [], options.imageSize, visFilters, []);
              combinedImg = uint8((modalImg + img) / 2);
              imwrite(modalImg, [outputImgFolder '/' fileName '_layer' num2str(vocabLevelItr) '_modal.png']);
              imwrite(combinedImg, [outputImgFolder '/' fileName '_layer' num2str(vocabLevelItr) '_combined.png']);
         end
         
         % Finally, save data structures.
         allNodes{vocabLevelItr} = nodes;
         allPrecisePositions{vocabLevelItr} = precisePositions;
         allActivations{vocabLevelItr} = nodeActivations;
    end
    
    numberOfInstances = sum(cellfun(@(x) numel(x), allActivations));
    
    %% If no instances have been found, exit.
    if numberOfInstances<1
        exportArr = [];
        return;
    end
    
    warning(w);
    
    %% Save all nodes to an array, and exit.
    activationArr = cat(1, allActivations{:});
    exportArr = zeros(numberOfInstances, 5, 'int32');
    exportArr(:,5) = 1;
    startIdx = 1;
    for vocabLevelItr = 1:numel(allNodes)
        thisLevelNodes = allNodes{vocabLevelItr};
        thisLevelPrecisePositions = allPrecisePositions{vocabLevelItr};
        if isempty(thisLevelNodes)
            break;
        end
        exportArr(startIdx:(startIdx + size(thisLevelNodes,1) - 1), 4) = vocabLevelItr;
        exportArr(startIdx:(startIdx + size(thisLevelNodes,1) - 1), 1) = thisLevelNodes(:,1);
        exportArr(startIdx:(startIdx + size(thisLevelNodes,1) - 1), 2:3) = thisLevelPrecisePositions;
        startIdx = startIdx + size(thisLevelNodes,1);
    end
end

