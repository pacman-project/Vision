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
function [exportArr, activationArr] = inferSubs(fileName, img, vocabulary, vocabularyDistributions, uniqueVocabularyChildren, vocabularyChildren, nodes, nodeActivations, options)
    % Read data into helper data structures.
    missingPartLog = log(0.00001);
    maxShareability = options.maxShareability;
    missingPartAllowed = false;
    poolFlag = true;
    labeledPooling = options.labeledPoolingTest;
    maxPosProb = 0.1;
    halfMatrixSize = (options.receptiveFieldSize+1)/2;
    maxSize = options.subdue.maxSize;
    if options.testDebug
         load([pwd '/filters/optimizationFilters.mat'], 'visFilters');
         outputImgFolder  = [options.currentFolder '/output/' options.datasetName '/reconstruction/test/' fileName];
         filterIds = round(((180/numel(options.filters)) * (0:(numel(options.filters)-1))) / (180/size(visFilters,3)))' + 1;
    end

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
    w = warning('off', 'all');
    
    %% Start discovery process.
    for vocabLevelItr = 2:min(numel(vocabulary), 10)    
         descriptions = double(uniqueVocabularyChildren{vocabLevelItr});
         prevLabelCounts = max(unique(prevVocabORLabels));
         multiplyVect = (double(prevLabelCounts+1).^((size(descriptions,2)-1):-1:0))';
         
         rfSize = getRFSize(options, vocabLevelItr);
         halfRealRFSize = floor(rfSize(1)/2) + 1;
         realRFSize = rfSize(1);
         rfRadius = options.receptiveFieldSizes(vocabLevelItr-1)-1;
         
         % Obtain edge novelty threshold.
         edgeNoveltyThr = max(options.minEdgeNoveltyThr, options.edgeNoveltyThr - options.edgeNoveltyThrRate * max(0, (vocabLevelItr-3)));
         edgeShareabilityThr = 1 - edgeNoveltyThr;
         if edgeShareabilityThr < 1
               allowedSharedLeafNodes = cellfun(@(x) numel(x), leafNodes) * edgeShareabilityThr;
         end
     
        % Calculate pool factor.
        poolFactor = nnz(~ismembc(2:vocabLevelItr, options.noPoolingLayers));
         
        %% Match subs from vocabLevel to their instance in graphLevel.
        % Obtain relevant data structures for the level.
        vocabLevel = vocabulary{vocabLevelItr};
        vocabLevelChildren = vocabularyChildren{vocabLevelItr};
        maxCount = size(vocabLevelChildren,2);
        
        vocabLevelDistributions = vocabularyDistributions{vocabLevelItr};
        if isempty(vocabLevel)
             break;
        end
        
        % Allocate space for this layer's realizations.
        layerNodes = cell(numel(vocabLevel),1);
        layerActivations = layerNodes;
        layerLeafNodes  = layerNodes;
        layerPrecisePositions = layerNodes;
        layerCombinations = layerNodes;
        layerChildren = layerNodes;
        
        % Obtain OR node labels.
        nodeORLabels = prevVocabORLabels(nodes(:,1));
        
        % Calculate pairwise distances.
        if size(nodes,1) == 1
             distances = false;
        else
             distances = squareform(pdist(single(nodes(:,2:3)), 'chebychev'));
             [X, Y]= meshgrid(1:size(nodes,1), 1:size(nodes,1));
             X = nodeORLabels(X);
             Y = nodeORLabels(Y);
             invalidEdges = X > Y;
             distances(1:(size(distances,1)+1):size(distances,1)^2) = Inf;
             distances = distances > 0 & distances <= rfRadius;
             distances(invalidEdges) = 0;
        end
        
        %% Apply leaf shareability filter and find single subs.
        loneNodeArr = ones(size(nodes,1),1) > 0;
%        loneNodeArr = zeros(size(nodes,1),1) > 0;
        adjMatrix = zeros(size(nodes,1)) > 0;
        for nodeItr = 1:size(nodes,1)
             % Check distance.
             adjacentNodes = find(distances(:, nodeItr));
             
             % Check edge matrices for filtering.
              if ~isempty(adjacentNodes)
                   if edgeShareabilityThr < 1
                        commonLeafCounts = cellfun(@(x) sum(ismembc(x, leafNodes{nodeItr})), leafNodes(adjacentNodes));
                        novelNodes = (commonLeafCounts <= allowedSharedLeafNodes(adjacentNodes)) & ...
                            (commonLeafCounts <= allowedSharedLeafNodes(nodeItr));
                        adjacentNodes = adjacentNodes(novelNodes);
                   end
              else
                   loneNodeArr(nodeItr) = true;
                   continue;
              end
             
              % If all adjacent nodes eliminated, move on.
              if isempty(adjacentNodes)
                   loneNodeArr(nodeItr) = true;
                   continue;
              else
                   adjMatrix(adjacentNodes, nodeItr) = 1;
              end
        end
        
        %% Find all node combinations.
        prevCliques = int32((1:size(nodes,1))');
        prevCliqueLabels = double(nodeORLabels(prevCliques));
        if size(prevCliques,1) == 1
             prevCliqueLabels = prevCliqueLabels';
        end
        
        allCliques = [];
        for sizeItr = 2:maxCount
               relevantDescriptions = descriptions(:, 1:sizeItr);
               relevantMultiplyVect = multiplyVect(((end+1)-sizeItr):end);
               relevantDescriptions = fastsortedunique((relevantDescriptions) * relevantMultiplyVect);
               prevCliqueLabels = prevCliqueLabels * relevantMultiplyVect(1:(end-1));
               
               nodeCount = size(adjMatrix,1);
               dummyLabels = (1:nodeCount)';

               % Find cliques of nodes. 
               % Allocate space for new cliques.
               newCliques = cell(size(prevCliques,1),1);

               % Find nodes adjacent to all nodes in a clique.
               for newCliqueItr = 1:size(prevCliques,1)
                    relevantRow = prevCliques(newCliqueItr, :);
                    newAdjNodes = int32(dummyLabels(all(adjMatrix(:, relevantRow), 2)));
                    
                    % Check to make sure only promising cliques are
                    % expanded.
                    newAdjNodes = newAdjNodes(ismembc(prevCliqueLabels(newCliqueItr) + double(nodeORLabels(newAdjNodes)), relevantDescriptions));
                    
                    % Save cliques.
                    if ~isempty(newAdjNodes)
                         newCliques{newCliqueItr} = cat(2, relevantRow(ones(size(newAdjNodes,1),1), :), newAdjNodes);
                    end
               end

               % Pad and save new cliques.
               newCliques = cat(1, newCliques{:});
               
               if isempty(newCliques)
                    break;
               end
               
               prevCliques = newCliques;
               prevCliqueLabels = double(nodeORLabels(prevCliques));
               if size(prevCliques,1) == 1
                    prevCliqueLabels = prevCliqueLabels';
               end
               newCliques = cat(2, newCliques, zeros(size(newCliques,1), maxSize - sizeItr, 1, 'int32'));
               allCliques = cat(1, allCliques, newCliques);
        end
%         
        %% Put statistics in cell arrays for faster reference.
         posArr = {vocabLevelDistributions.childrenPosDistributionProbs};
         minPosLogProbs = cat(1, vocabLevelDistributions.minPosActivationLog);
         minLogProbs = cat(1, vocabLevel.minActivationLog);
         childrenCounts = sum(vocabLevelChildren>0,2);
         
        %% Starting inference from layer 2. 
         if ~isempty(allCliques)
              labeledCliques = allCliques;
              labeledCliques(labeledCliques > 0) = nodeORLabels(labeledCliques(labeledCliques > 0));
              labeledCliques = labeledCliques(:, 1:maxCount);

              % Eliminate some cliques.
              labeledCliques = double(labeledCliques) * multiplyVect;
              vocabLevelChildren = double(vocabLevelChildren) * multiplyVect;

              for vocabItr = 1:numel(vocabLevel)
                    % Find cliques.
                    instanceCombinations = allCliques(labeledCliques == vocabLevelChildren(vocabItr), 1:childrenCounts(vocabItr));

                    % Single node elimination.
     %                if numel(vocabNodeChildren) == 1
     %                     validInstances = loneNodeArr & validInstances;
     %                end


                    % If no centers have matched, continue.
                    if isempty(instanceCombinations)
                      continue;
                    end

                    numberOfCombinations = size(instanceCombinations,1);
                    numberOfChildren = size(instanceCombinations,2);

                    if numberOfCombinations == 0
                         continue;
                    end

                   %% Check leaf node shareability, and save leaf nodes for every instance.
                   instanceLeafNodes = cell(numberOfCombinations, 1);
                   if edgeShareabilityThr < 1
                        for instanceItr = 1:numberOfCombinations
                             tempCombinations = instanceCombinations(instanceItr,:);
                             tempArr = sort(cat(1, leafNodes{tempCombinations(~isinf(tempCombinations))}));
                             tempLeafNodes = fastsortedunique(tempArr);
                             instanceLeafNodes{instanceItr} = tempLeafNodes;
                        end
                   end

                   %% Combinations are discovered. Now, we calculate an activation value for every combination.
                   % First, we have to pick relevant distributions for each node. We
                   % need combined relative positions first.
                   % Calculate instance positions.
                   instancePrecisePositions = zeros(numberOfCombinations,2, 'single');
                   for instanceItr = 1:numberOfCombinations
                       positions = precisePositions(instanceCombinations(instanceItr, :)',:);
                       instancePrecisePositions(instanceItr,:) = round((min(positions) + max(positions))/2);
                   end

                  %% Calculate activations!
                  childrenProbs = posArr{vocabItr};
                  childrenPosActivations = zeros(size(instanceCombinations), 'single');
                  for childItr = 1:size(instanceCombinations,2)
                       relevantChildren = instanceCombinations(:, childItr);
                       relevantProbs = childrenProbs{childItr};

                       % Calculate position probabilities.
                       samples = (precisePositions(relevantChildren, :) - instancePrecisePositions) + halfRealRFSize;
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
                  if childrenCounts(vocabItr) == 1
                        layerChildren{vocabItr} = instanceCombinations(validCombinations,:);
                  else
                        layerCombinations{vocabItr} = instanceCombinations(validCombinations,:);
                  end
              end
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
               if labeledPooling
                    [~, IA, ~] = unique(combinedArr(:,2:4), 'rows', 'stable');
               else
                    [~, IA, ~] = unique(combinedArr(:,3:4), 'rows', 'stable');
               end
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

