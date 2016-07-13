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
    missingPartAllowed = false;
    poolFlag = true;
 %   activationThr = log(0.001);
    activationThr = -10;
    dummyPosSigma = 0.001;
    maxPosProb = 0.1;
    halfMatrixSize = (options.receptiveFieldSize+1)/2;
    maxLevel = 10;
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
    allActivations(1) = {nodeActivations};
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
     
    for vocabLevelItr = 2:min(numel(vocabulary), 10)    
         realRFSize = getRFSize(options, vocabLevelItr);
         halfRealRFSize = round(realRFSize(1)/2);  
        
         modes = allModes{vocabLevelItr-1};
         curModeProbs = modeProbs{vocabLevelItr-1};
         rankModes = int32(unique(modes(:,1:2), 'rows'));
         
         if ismember(vocabLevelItr-1, options.smallRFLayers)
             rfRadius = smallHalfMatrixSize;
         else
             rfRadius = halfRFSize;
         end
         
        % If no pooling has been performed at this layer, and previous layer
        % has small RF, we have a minimum RF limit. 
        if ismember(vocabLevelItr-1, options.noPoolingLayers) && ismember(vocabLevelItr - 2, options.smallRFLayers) && vocabLevelItr < 8
            minRFRadius = floor(halfMatrixSize/2);
        else
            minRFRadius = 1;
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
        
        %% Starting inference from layer 2. 
         for vocabItr = 1:numel(vocabLevel)
             % Match the center first.
             vocabNode = vocabLevel(vocabItr);
             vocabNodeChildren = vocabNode.children;
             vocabNodeDistributions = vocabLevelDistributions(vocabItr);
             centerId = vocabNodeChildren(1);
             validInstances = nodeORLabels == centerId;
             
             % Obtain activation thresholds.
             minActivationLog = vocabNode.minActivationLog;
             minPosActivationLog = vocabNodeDistributions.minPosActivationLog;
             
             if numel(vocabNodeChildren) == 1
                  continue;
             end
             
             % Obtain edge matrices.
             if numel(vocabNodeChildren) > 1
                   secChildren = vocabNodeChildren(2:end)';
                   matrixIdx = zeros(numel(secChildren),1);
                   for itr = 1:numel(matrixIdx)
                      matrixIdx(itr) = find(rankModes(:,1) == vocabNodeChildren(1) & rankModes(:,2) == secChildren(itr));
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
             instanceCombinations = cell(numberOfCenterNodes,1);

               %% Iteratively match instances and get all node combinations.
               for centerItr = 1:numberOfCenterNodes
                       instanceChildren(centerItr, 1) = {centerNodes(centerItr)};

                       validCols = ones(numel(vocabNodeChildren),1) > 0;
                       for childItr = 2:numel(vocabNodeChildren)
                             tempNodes = find(nodeORLabels == vocabNodeChildren(childItr) & distances(:, centerNodes(centerItr)) < rfRadius & distances(:, centerNodes(centerItr)) >= minRFRadius);
                             
                             if ~isempty(tempNodes)
                                  % Fine-tune the system by filtering temp nodes
                                  % more.
                                  curEdgeType = vocabNode.adjInfo(childItr-1, 3);
                                  curModes = modes(modes(:,1) == vocabNodeChildren(1) & modes(:,2) == vocabNodeChildren(childItr), :);
                                  relativeCoords = nodes(tempNodes,2:3) - repmat(nodes(centerNodes(centerItr), 2:3), numel(tempNodes),1) + halfRFSize;
                                  linIdx = relativeCoords(:,1) + (relativeCoords(:,2)-1)*RFSize;
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
               
              %% Save leaf nodes.
              instanceLeafNodes = cell(numberOfCombinations, 1);
              for instanceItr = 1:numberOfCombinations
                   tempCombinations = instanceCombinations(instanceItr,:);
                   tempLeafNodes = fastsortedunique(sort(cat(1, leafNodes{tempCombinations(~isinf(tempCombinations))})));
                   instanceLeafNodes{instanceItr} = tempLeafNodes;
              end
               
              %% Combinations are discovered. Now, we calculate an activation value for every combination.
              % First, we have to pick relevant distributions for each node. We
              % need combined relative positions first.
              
              % Calculate instance positions.
              instancePrecisePositions = zeros(numberOfCombinations,2, 'single');
              for instanceItr = 1:numberOfCombinations
                  instancePrecisePositions(instanceItr,:) = round(mean(leafPrecisePositions(instanceLeafNodes{instanceItr},:), 1));
              end
             
             %% Fill in joint positions.
             if numberOfChildren > 1
%                  instanceCombinationLabels = instanceCombinations;
%                   validChildIdx = ~isinf(instanceCombinationLabels);
%                  instanceCombinationLabels(validChildIdx) = nodes(instanceCombinationLabels(validChildIdx),1);
                  
                  instanceChildrenCombinedPos = zeros(size(instanceCombinations,1), size(instanceCombinations,2) * 2);
                  % Put children positions in an array.
                  for instanceItr = 1:size(instanceCombinations,2)
                       relevantChildren = instanceCombinations(:, instanceItr);

                       % If we're working with peripheral sub-parts, we calculate
                       % position distributions as well.
                       samples = precisePositions(relevantChildren, :) - instancePrecisePositions;

                       % Save samples (positions and labels).
                       instanceChildrenCombinedPos(:, ((instanceItr-1)*2+1):((instanceItr)*2)) = samples;
                  end
%                  [~, instanceJointPos] = vocabNodeDistributions.predictMissingInfo(instanceCombinationLabels, instanceCombinations, precisePositions);
                  validCombinations = ones(numberOfCombinations,1) > 0;
                  
                  %% Calculate activations!
                  % Calculate position likelihood.
                  % TODO: Encapsulate with try/catch.
%                  try
                   posProbs = pdf(vocabNodeDistributions.childrenPosDistributions, instanceChildrenCombinedPos / halfRealRFSize); 
                   activationDenom = max(pdf(vocabNodeDistributions.childrenPosDistributions, vocabNodeDistributions.childrenPosDistributions.mu)) * 1/maxPosProb;
                   posActivations = posProbs / activationDenom;
                   posActivations = single(log(posActivations));
                   
                   % Apply a position likelihood-based inhibition here.
                   validCombinations = validCombinations & posActivations >= minPosActivationLog;

                   % Calculate overall activations.
                   prevChildrenActivations = nodeActivations(instanceCombinations);
                   if ~isequal(size(prevChildrenActivations), size(instanceCombinations))
                       prevChildrenActivations = prevChildrenActivations';
                   end
                   avgActivations = mean([prevChildrenActivations, posActivations], 2);
                   
                   % Apply overall activation threshold here.
                   validCombinations = validCombinations & avgActivations >= minActivationLog;
                   instanceActivations = avgActivations;
                  
             else
                  instanceActivations = nodeActivations(instanceCombinations);
                  validCombinations = instanceActivations >= minActivationLog; 
             end

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
         experts = projectNode([nodes(:,1), precisePositions, ones(size(nodes,1),1)*vocabLevelItr], vocabularyDistributions, 'modal', options);
         experts(:,1) = filterIds(experts(:,1));
         modalImg = obtainPoE(experts, [], [], options.imageSize, visFilters, []);
         combinedImg = uint8((modalImg + img) / 2);
  %       figure, imshow(combinedImg);
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

