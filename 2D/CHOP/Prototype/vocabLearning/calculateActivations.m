%> Name: calculateActivations
%>
%> Description: Calculates activation values for the data and given
%> realizations. For each realization, the level 1 parts which are not
%> explained are penalized with a small probability.
%> 
%> @param graphLevel List of realizations to calculate activations of. 
%> @param level1Coords Level 1 nodes, with reduced coordinates. 
%>
%> @retval extendedSubs Extended sub list.
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 18.01.2016
function [ graphLevel ] = calculateActivations( vocabLevel, vocabulary, graphLevel, mainGraph, level1Coords, options, levelItr )
     % Program variables.
     epsilon = 0.0001;
     logEpsilon = single(log(epsilon));
     poolDim = options.poolDim;  
     if strcmp(options.filterType, 'gabor')
          stride = options.gabor.stride;
          inhibitionRadius = options.gabor.inhibitionRadius;
     else
          stride = options.auto.stride;
          inhibitionRadius = options.auto.inhibitionRadius;
     end

     
     % First, we put level 1 nodes into bins, each of which belongs to an
     % image.
     numberOfImages = max(level1Coords(:,1));
     imageOffsets = zeros(numberOfImages+1,1);
     for imageItr = 1:numberOfImages
          idx = level1Coords(:,1) == imageItr;
          startItr = find(idx,1, 'first');
          imageOffsets(imageItr) = startItr;
     end
     
     % First, we create some intermediate data structures for speed.
     mainGraph{levelItr} = graphLevel;
     allPositions = cell(levelItr,1);
     allPrecisePositions = cell(levelItr,1);
     allRealLabelIds = cell(levelItr,1);
     for curLevelItr = 1:levelItr
        allPositions{curLevelItr} = single(cat(1, mainGraph{curLevelItr}.position));
        allPrecisePositions{curLevelItr} = single(cat(1, mainGraph{curLevelItr}.precisePosition));
        allRealLabelIds{curLevelItr} = single([mainGraph{curLevelItr}.realLabelId]);
     end
     
    % Obtain leaf nodes that (theoretically) could be covered by every
    % instance.
    prevLevel = mainGraph{levelItr-1};
    prevPositions = cat(1, prevLevel.position);
    possibleLeafNodes = cell(size(prevPositions,1),1);
    bounds = calculateRFBounds(prevPositions, levelItr-1, options, false);
    for nodeItr = 1:numel(prevLevel)
        possibleLeafNodes{nodeItr} = int32(find(level1Coords(:,1) == prevLevel(nodeItr).imageId & ...
             level1Coords(:,2) >= bounds(nodeItr,1) & level1Coords(:,2) <= bounds(nodeItr,3) & ...
             level1Coords(:,3) >= bounds(nodeItr,2) & level1Coords(:,3) <= bounds(nodeItr,4)));
    end
     
     % Go through the list of realizations, and find the activation value
     % for each of them. 
     parfor nodeItr = 1:numel(graphLevel)
          logLikelihood = 0;
          decisionCtr = 0;
          curLevelItr = levelItr;
          node = graphLevel(nodeItr);
          
          % Add penalty for missing data.
          missingNodeCount = numel(possibleLeafNodes{node.children(1)}) - numel(node.leafNodes);
          logLikelihood = logLikelihood + missingNodeCount * logEpsilon;
          decisionCtr = decisionCtr + missingNodeCount;
          
          % Calculate likelihood of the tree.
          % For each vocabulary node, we generate lots and lots of samples.
          % Then, we perform kernel density estimation (simply put, how
          % many samples fell into the defined area?). That gives us the
          % likelihood of the data given the mode. It is approximate, but
          % the approximation should get better and better as we generate
          % more samples. In 2D case, we can calculate exact probability by
          % integral. 
          nodes = node;
          
         %% First, we recursively calculate probabilities of the previous layer's choices. 
         curVocabLevel = vocabLevel;
         while curLevelItr > 1.001
             newNodes = cell(size(nodes,1),1);
             prevGraphLevel = mainGraph{curLevelItr-1};
%             positions = allPositions{curLevelItr-1};
             precisePositions = allPrecisePositions{curLevelItr-1};
             flexRadius = floor(options.receptiveFieldSize/2)^(curLevelItr-1) * stride * (inhibitionRadius+1);
             realLabelIds = allRealLabelIds{curLevelItr-1};
             for itr = 1:size(nodes,1)
                 vocabNode = curVocabLevel(nodes(itr).realLabelId);

                 % Sample from the discrete and continuous distributions.
                 childrenLabelDistributions = vocabNode.childrenLabelDistributions;
                 labelSamples = realLabelIds(nodes(itr).children);
                 probIdx = ismember(childrenLabelDistributions(:,1:(end-1)), labelSamples, 'rows');
                 if nnz(probIdx) == 1
                      labelProb = childrenLabelDistributions(probIdx, end);
                      decisionCtr = decisionCtr + 1;
                 else
                      labelProb = logEpsilon;
                      decisionCtr = decisionCtr + 1;
                      logLikelihood = logLikelihood + labelProb;
                      continue;
                 end

                 % We calculate the probability of a sample given the node
                 % distributions. If D = 2 (One sub-child), we simply
                 % calculate integral over CDF. In case of more
                 % sub-children, we simply use a density estimator.
                 posDistributions = vocabNode.childrenPosDistributions{1};
                 posDistributions = posDistributions{probIdx};
                 posSamples = vocabNode.childrenPosSamples{1};
                 posSamples = posSamples{probIdx};
                 if ~isempty(posDistributions)
                      sampledPositions = precisePositions(nodes(itr).children, :) - single(repmat(precisePositions(nodes(itr).children(1), :), numel(nodes(itr).children), 1));
                      
                      % From these coarse positions, let's find boundaries
                      % of the RFs in precise position space. We add 1 to
                      % sample positions to account for the fact that in
                      % images, indexing starts from 1, not 0.
%                      sampledPositions = calculateRFBounds(sampledPositions+1, levelItr-1, options, true) -1;
 %                     sampledPositions = sampledPositions - repmat(meanPos, 1, 2);
                      sampledPositions = sampledPositions(2:end,:);
                      sampledPositions = sampledPositions';
                      sampledPositions = sampledPositions(:)';
                      
%                       sampledPrecisePositions = sampledPrecisePositions(2:end,:);
%                       sampledPrecisePositions = sampledPrecisePositions';
%                       sampledPrecisePositions = sampledPrecisePositions(:)';
                      
                      % From these sampled positions, we get back to the
                      % precise position space and select a boundary for
                      % calculating probabilities.
                      % First, we take into account the pooling so far.
                      sampledPositionsStart = sampledPositions - flexRadius;
                      sampledPositionsEnd = sampledPositions + flexRadius;
                      
                      % Sanity check.
%                       checkIdx = sum(sampledPrecisePositions >= repmat(sampledPositionsStart, size(sampledPrecisePositions,1),1) & ...
%                           sampledPrecisePositions < repmat(sampledPositionsEnd, size(sampledPrecisePositions,1),1),2) == size(posSamples,2);
%                       if nnz(checkIdx) ~= numel(checkIdx)
%                            error('Error in calculateActivations: Check for positions failed.');
%                       end
                      
                      % Multiple sub-parts.
                      sampleCount = size(posSamples,1);
                      validIdx = sum(posSamples >= repmat(sampledPositionsStart, sampleCount,1) & ...
                           posSamples < repmat(sampledPositionsEnd, sampleCount,1),2) == size(posSamples,2);
                      posProb = nnz(validIdx) / sampleCount;
                 else
                      posProb = single(1);
                 end
                 decisionCtr = decisionCtr + 1;
                 
                 % Update log likelihood.
                 logLikelihood = log(labelProb) + max(logEpsilon, log(posProb)) + logLikelihood;

                 % Finally, obtain new nodes and calculate likelihood for
                 % them.
                 newNodeSet = prevGraphLevel(nodes(itr).children);
                 newNodes{itr} = newNodeSet;
             end
             nodes = cat(1, newNodes{:});
             curLevelItr = curLevelItr - 1;
             curVocabLevel = vocabulary{curLevelItr};
         end
         
         % Calculate activations.
         avgNegLogLikelihood = - (logLikelihood / decisionCtr);
         graphLevel(nodeItr).activation = (abs(logEpsilon) - avgNegLogLikelihood) / (abs(logEpsilon));
     end
end

