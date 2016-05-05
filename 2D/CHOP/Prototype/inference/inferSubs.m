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
    halfSize = ceil(RFSize/2);
    missingPartLog = log(0.00001);
    missingPartAllowed = false;
 %   activationThr = log(0.001);
    activationThr = -Inf;
    dummyPosSigma = 1;
    maxLevel = 8;
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
%    allLeafNodes = cell(numel(vocabulary),1);
    
    % Fill layer 1 data structures.
    precisePositions = single(nodes(:,4:5));
    leafNodes = num2cell((1:size(nodes,1))');
    nodeActivations = log(nodeActivations);
    
    allNodes(1) = {nodes(:,1:3)};
    allActivations(1) = {nodeActivations};
%    allLeafNodes{1} = leafNodes;
    allPrecisePositions{1} = precisePositions;
    prevVocabORLabels = cat(1, vocabulary{1}.label);
    
    % Get RF choices (to normalize position likelihood).
    rfChoices = RFSize^2;
     
    for vocabLevelItr = 2:min(numel(vocabulary), maxLevel)
         modes = allModes{vocabLevelItr-1};
         curModeProbs = modeProbs{vocabLevelItr-1};
         rankModes = unique(modes(:,1:2), 'rows');
         
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
        vocabNodeLabels = cat(1, vocabLevel.label);
        
        % Calculate pairwise distances.
        if size(nodes,1) == 1
             distances = Inf;
        else
             distances = squareform(pdist(single(nodes(:,2:3))));
             distances(1:(size(distances,1)+1):size(distances,1)^2) = Inf;
        end
        
        %% Starting inference from layer 2. 
         for vocabItr = 1:numel(vocabLevel)
             % Match the center first.
             vocabNode = vocabLevel(vocabItr);
             vocabNodeChildren = vocabNode.children;
             vocabNodeDistributions = vocabLevelDistributions(vocabItr);
             centerId = vocabLevel(vocabItr).children(1);
             validInstances = nodeORLabels == centerId;
             
             if numel(vocabNodeChildren) == 1
                  continue;
             end
             
             % Obtain edge matrices.
             if numel(vocabNodeChildren) > 1
                   secChildren = vocabNodeChildren(2:end)';
                   labelPairs = [ones(numel(secChildren),1, 'int32') * vocabNodeChildren(1), secChildren];
                   [~, matrixIdx] = ismember(labelPairs, rankModes, 'rows');
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
             instanceCombinations = cell(numberOfCenterNodes,1);

               %% Iteratively match instances and get all node combinations.
               for centerItr = 1:numberOfCenterNodes
                       instanceChildren(centerItr, 1) = {centerNodes(centerItr)};

                       validCols = ones(numel(vocabNode.children),1) > 0;
                       for childItr = 2:numel(vocabNode.children)
                             tempNodes = find(nodeORLabels == vocabNodeChildren(childItr) & distances(:, centerNodes(centerItr)) < halfSize);
                             
                             if ~isempty(tempNodes)
                                  % Fine-tune the system by filtering temp nodes
                                  % more.
                                  curEdgeType = vocabNode.adjInfo(childItr-1, 3);
                                  curModes = modes(modes(:,1) == vocabNodeChildren(1) & modes(:,2) == vocabNodeChildren(childItr), :);
                                  relativeCoords = nodes(tempNodes,2:3) - repmat(nodes(centerNodes(centerItr), 2:3), numel(tempNodes),1) + halfSize;
                                  linIdx = sub2ind([RFSize, RFSize], relativeCoords(:,1), relativeCoords(:,2));
                                  curEdgeMatrix = edgeMatrices(:,:,childItr-1);
                                  
                                  % Save edge types.
                                  edgeMatrixIndices = curEdgeMatrix(linIdx);
                                  edgeTypes = zeros(size(edgeMatrixIndices),'single');
                                  if find(edgeMatrixIndices) > 0
                                       edgeTypes(edgeMatrixIndices>0) = curModes(edgeMatrixIndices(edgeMatrixIndices > 0), 3);
                                  end
                                  tempNodes = tempNodes(edgeTypes == curEdgeType);
                             end
                             
                             % If tempNodes is not empty, save it.
                             if isempty(tempNodes)
                                  tempNodes = Inf;
                                  validCols(childItr) = 0;
                             end
                             instanceChildren(centerItr, childItr) = {tempNodes};
                       end

                       % Find all node combinations (not containing node
                       % repetitions).
                       curInstanceCombinations = allcomb(instanceChildren{centerItr,validCols});
                       sortedmatrix = sort(curInstanceCombinations,2);
                       validInstanceIdx = all(diff(sortedmatrix,[],2)~=0,2);
                       if nnz(validCols) ~= numel(validCols)
                            curInstanceCombinations = allcomb(instanceChildren{centerItr,:});
                       end
                       curInstanceCombinations = curInstanceCombinations(validInstanceIdx,:);
                       instanceCombinations{centerItr} = curInstanceCombinations;
               end

               % We have the instances for this sub. Save the info.
               instanceCombinations = cat(1, instanceCombinations{:});
               
               % Remove rows having inf, if prediction is not needed.
              if ~missingPartAllowed
                   instanceCombinations = instanceCombinations(~any(isinf(instanceCombinations),2),:);
              end
               
               numberOfCombinations = size(instanceCombinations,1);
               numberOfChildren = size(instanceCombinations,2);
               
               if numberOfCombinations == 0
                    continue;
               end
             %% Combinations are discovered. Now, we calculate an activation value for every combination.
             % First, we have to pick relevant distributions for each node. We
             % need combined relative positions first.
             % Learn mean positions of the children.
             instancePrecisePositions = precisePositions(instanceCombinations(:,1), :);
             
             %% Fill in joint positions.
             if numberOfChildren > 1
                  instanceCombinationLabels = instanceCombinations;
                  validChildIdx = ~isinf(instanceCombinationLabels);
                  instanceCombinationLabels(validChildIdx) = nodes(instanceCombinationLabels(validChildIdx),1);
                  [~, instanceJointPos] = vocabNodeDistributions.predictMissingInfo(instanceCombinationLabels, instanceCombinations, precisePositions);

                  %% Calculate activations!
                  % Calculate position likelihood.
                  % TODO: Encapsulate with try/catch.
                  try
                       instanceActivations = log(pdf(vocabNodeDistributions.childrenPosDistributions, instanceJointPos) / rfChoices);
                  catch %#ok<CTCH>
                       % Calculate pseudo-probability based on distance to
                       % the peaks.
                       tempDist = gmdistribution(vocabNodeDistributions.childrenPosDistributions.mu(1,:), ones(1, (numberOfChildren-1)*2)*dummyPosSigma);
                       instanceActivations = log(pdf(tempDist, instanceJointPos) / rfChoices);
                  end

                  % First, we need to calculate precise positions.
                  for instanceItr = 1:numberOfCombinations
                      centerCoords = precisePositions(instanceCombinations(instanceItr,1), :);
                      addedCoords = repmat(centerCoords, (numberOfChildren-1), 1);
                      secChildrenPos = reshape(instanceJointPos(instanceItr, :), 2, (numberOfChildren-1))';
                      secChildrenPos = secChildrenPos + addedCoords;
                      childrenPos = cat(1, centerCoords, secChildrenPos);
                      precisePosition = round((min(childrenPos,[], 1) + max(childrenPos, [], 1)) / 2);
                      instancePrecisePositions(instanceItr, :) = precisePosition;

                      % Get activations of true children, add them to the sum.
                      trueChildren = instanceCombinations(instanceItr, :);
                      trueChildren = trueChildren(~isinf(trueChildren));
                      numberOfMissingChildren = numberOfChildren - numel(trueChildren);
                      instanceActivations(instanceItr) = mean([instanceActivations(instanceItr); nodeActivations(trueChildren); ones(numberOfMissingChildren, 1) * missingPartLog]);
                  end
             else
                  instanceActivations = nodeActivations(instanceCombinations);
             end

             %% Save leaf nodes.
             instanceLeafNodes = cell(numberOfCombinations, 1);
             for instanceItr = 1:numberOfCombinations
                   tempCombinations = instanceCombinations(instanceItr,:);
                   tempLeafNodes = cat(1, leafNodes{tempCombinations(~isinf(tempCombinations))});
                   instanceLeafNodes{instanceItr} = tempLeafNodes;
             end

             %% Create new nodes.
             newNodes = [ones(numberOfCombinations, 1, 'int32') * vocabItr, int32(instancePrecisePositions)];

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
         
         if isempty(nodes)
              break;
         end
         
         %% Perform likelihood-based inhibition.
         validActivations = nodeActivations>= activationThr;
         nodes = nodes(validActivations,:);
         nodeActivations = nodeActivations(validActivations,:);
         leafNodes = leafNodes(validActivations, :);
         precisePositions = precisePositions(validActivations, :);
         pooledPositions = calculatePooledPositions(precisePositions, vocabLevelItr+1, options);
         
         if isempty(nodes)
              break;
         end
         
          % Get unique nodes for each label, coords, activation triplet.
          combinedArr = double([double(vocabNodeLabels(nodes(:,1))), pooledPositions, -nodeActivations]);

          % First, we order the nodes.
          [combinedArr, sortIdx] = sortrows(combinedArr);
          
          %% Perform max pooling and save data structures.
          [~, IA, ~] = unique(combinedArr(:,1:3), 'rows', 'stable');

          % Save remaining data and move on to the next layer.
         validIdx = sort(sortIdx(IA));
         nodes = nodes(validIdx,:);
         nodeActivations = nodeActivations(validIdx,:);
         leafNodes = leafNodes(validIdx, :);
         precisePositions = precisePositions(validIdx, :);
         prevVocabORLabels = cat(1, vocabLevel.label);
         nodes(:,2:3) = calculatePooledPositions(precisePositions, vocabLevelItr+1, options);
         
         %% For debugging purposes, we visualize the nodes.
         experts = projectNode([nodes(:,1), precisePositions, ones(size(nodes,1),1)*vocabLevelItr], vocabularyDistributions, 'modal');
         experts(:,1) = filterIds(experts(:,1));
         modalImg = obtainPoE(experts, [], [], options.imageSize, visFilters, []);
         combinedImg = uint8((modalImg + img) / 2);
         figure, imshow(combinedImg);
         imwrite(modalImg, [outputImgFolder '/' fileName '_layer' num2str(vocabLevelItr) '_modal.png']);
         imwrite(combinedImg, [outputImgFolder '/' fileName '_layer' num2str(vocabLevelItr) '_combined.png']);
         
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

