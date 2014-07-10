%> Name: combineParts
%>
%> Description: Combines parts in a single layer if they essentially encode
%> the same representation. This function helps CHOP to
%> generalize over the dataset better. Parts are less precise, but more
%> general.
%>
%> @param vocabLevel The current vocabulary level, ASSUMED ordered by
%> mdlValues.
%> @param graphLevel The current graph level, ASSUMED ordered by labelIds.
%> @param currentModes Current modes representing geometric relationships.
%> @param similarityMatrix Part similarity matrix of the previous layer.
%> @param options Program options.
%> 
%> @retval vocabLevel Modified vocab level.
%> @retval graphLevel Modified graph level.
%> @retval newSimilarityMatrix New similarity matrix of the next level.
%> @retval subClasses Array containing new label id of each composition. 
%> Some compositions retain their ids, while others take ids of persistent
%> compositions.
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 05.05.2014
function [vocabLevel, graphLevel, newSimilarityMatrix, subClasses] = combineParts(vocabLevel, graphLevel, currentModes, similarityMatrix, options)
    global simMat;
    simMat = similarityMatrix;
    regularizationParam = options.subdue.maxSize * 2;
    threshold = options.subdue.threshold * regularizationParam;
    coarsityParam = 4;
    %% First inhibition is done via structural evaluation. 
    % If the structures of two different subs match, it is very likely that
    % they have same node sets. We must eliminate such cases.
    % We only need mode positions.
    if options.useReceptiveField
        %% If no current modes exist, we eliminate them based on their children only.
        if isempty(currentModes)
            vocabDescriptions = {vocabLevel.children};
            vocabDescriptions = cellfun(@(x) mat2str(sort(x)), vocabDescriptions, 'UniformOutput', false);
            [~, IA, IC] = unique(vocabDescriptions, 'stable');
            vocabLevel = vocabLevel(1, IA);
            
            % Update graph level.
            labelIds = cat(1, graphLevel.labelId);
            IC = num2cell(IC(labelIds));
            [graphLevel.labelId] = deal(IC{:});
        else
            %% If current modes do exist, we can use the positioning information to eliminate parts further.
            currentModes = currentModes(:,3:4);
            currentModes((size(currentModes,1)+1),:) = [0, 0];

            % Get unique sub descriptions independent of the center node. Each
            % description is a string that encodes coarse structure of the sub.
            vocabNodeLabels = {vocabLevel.children};
            vocabEdges = {vocabLevel.adjInfo};
            newMode = size(currentModes,1);
            vocabNeighborModes = cellfun(@(x) [newMode; x(:,3)], vocabEdges, 'UniformOutput', false);
            vocabNodePositions = cellfun(@(x) currentModes(x,:) - repmat(min(currentModes(x,:)), numel(x), 1), vocabNeighborModes, 'UniformOutput', false);

            % Get rid of noise in node positions to have a coarse geometric
            % representation.
%            normalizeConstant = coarsityParam/max(sqrt(sum(currentModes.^2,2)));
            normalizeConstant = 1;
            
            % Maximum distance between two nodes in 2D space.
            maxDistance = max(sqrt(sum(currentModes.^2,2)))*2;
            
            % Downsample node positions so we have coarse representations.
            vocabNodePositions = cellfun(@(x) round(x*normalizeConstant), vocabNodePositions, 'UniformOutput', false);

            vocabSortOrder = cell(size(vocabLevel,1),1);
            for vocabNodeItr = 1:numel(vocabLevel)
                [~, vocabSortOrder{vocabNodeItr}] = sortrows(vocabNodePositions{vocabNodeItr});
            end
            vocabDescriptions = cellfun(@(x,y,z) [x(z)', y(z,:)], vocabNodeLabels, vocabNodePositions, vocabSortOrder, 'UniformOutput', false);
%            vocabDescriptions = cellfun(@(x) mat2str(x(:)'), vocabDescriptions, 'UniformOutput', false);

            % Now, we have unique vocabulary descriptions, independent of center
            % nodes. Eliminate those match, by keeping only one (the one with
            % highest mdlScore) of each set.
     %       vocabDescriptions = vocabDescriptions(numel(vocabLevel):-1:1);
%            [~, IA, IC] = unique(vocabDescriptions, 'stable');
            
            %% For each sub, we will find matching ones with a cost.
            numberOfSubs = numel(vocabDescriptions);
            matchedSubs = zeros(numberOfSubs,1);
            subClasses = (1:numberOfSubs)';
            newSimilarityMatrix = zeros(numberOfSubs);
            
            % Here, we parallelize the part comparison process. In order to
            % do that, numberOfThreads best parts are selected, while the
            % rest of the parts are compared with the selected ones. When
            % comparison is over, the data structures are updated, and a
            % new set of parts are selected from remaining list.
 %           numberOfThreads = options.numberOfThreads;
%             partStartIdx = 1;
%  %           newSimilarityMatrix = newSimilarityMatrix(:);
%             
%             while true
%                 % Select first numberOfThreads parts.
%                 [selectedParts, partEndIdx] = GetBestParts(vocabDescriptions, partStartIdx, ...
%                     numberOfThreads, maxDistance, threshold);
%                 similarityMatrixEntries = cell(numberOfSubs-partStartIdx,1);
%                 
%                 for partItr = (partStartIdx+1):numberOfSubs
%                     description = vocabDescriptions{partItr};
%                     
%                     newEntries = zeros(1, numel(selectedParts));
%                     for partItr2 = 1:numel(selectedParts)
%                         description2 = vocabDescriptions{selectedParts(partItr2)}; %#ok<PFBNS>
%                         matchingCost = InexactMatch(description, description2, maxDistance);
%                         if matchingCost <= threshold && selectedParts(partItr2) < partItr
%                             matchedSubs(partItr) = 1;
%                             subClasses(partItr) = selectedParts(partItr2);
%                             break; %#ok<PFBR>
%                         else
%                             newEntries(partItr2) = matchingCost;
%                         end
%                     end
%                     if ~matchedSubs(partItr)
%                         similarityMatrixEntries{partItr} = newEntries;
%                     end
%                 end
%                 
%                 % Update the similarity matrix based on the new entries.
%                 allEntries = cat(1, similarityMatrixEntries{:});
%                 newSimilarityMatrix((partStartIdx+1):numberOfSubs, selectedParts) = allEntries;
%                 newSimilarityMatrix(selectedParts, (partStartIdx+1):numberOfSubs) = allEntries';
%                 
%                 partStartIdx = partEndIdx;
%             end
            
            for subItr = 1:(numberOfSubs-1)
                if ~matchedSubs(subItr) 
                    % Get the description of this composition.
                    description = vocabDescriptions{subItr};
                    
                    % Run inexact matching on every pair of subs.
                    for subItr2 = (subItr+1):numberOfSubs
                        if ~matchedSubs(subItr2)
                            description2 = vocabDescriptions{subItr2};
                            
                            % Match the two descriptions in an inexact
                            % manner. TODO: This part will be improved.
                            matchingCost = InexactMatch(description, description2, maxDistance);
                            
                            % If they match, assign the class of the
                            % first composition to that of the matching
                            % one. TODO: This part will be improved to
                            % assign each composition to its closest match,
                            % not highest scoring match.
                            if matchingCost <= threshold
                                matchedSubs(subItr2) = 1;
                                subClasses(subItr2) = subItr;
                            else
                                newSimilarityMatrix(subItr, subItr2) = matchingCost;
                                newSimilarityMatrix(subItr2, subItr) = matchingCost;
                            end
                        end
                    end
                end
            end
            
            % Process the last sub, if it has not been assigned to any
            % other sub.
            if ~matchedSubs(end)
                subClasses(end) = numberOfSubs;
            end
            
            % Create new similarity matrix.
%            newSimilarityMatrix = newSimilarityMatrix(remainingClasses, remainingClasses);
            newSimilarityMatrix = newSimilarityMatrix/max(max(newSimilarityMatrix));
            
            % Update graph level's ids, and sort them based on labelIds.
            labelIds = subClasses(cat(1, graphLevel.labelId));
            imageIds = cat(1, graphLevel.imageId);
            arrayToSort = [imageIds, labelIds];
            [~, sortedIdx] = sortrows(arrayToSort);
            IC = num2cell(labelIds);
            [graphLevel.labelId] = deal(IC{:});
            graphLevel = graphLevel(sortedIdx);
        end
    end
end

%> Name: GetBestParts
%>
%> Description: 
%>
%> @param 
%>
%> @retval 
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 10.07.2014
function [selectedParts, partEndIdx] = GetBestParts(vocabDescriptions, partStartIdx, ...
                    numberOfThreads, maxDistance, threshold)
    selectedPartCount = 1;
    partEndIdx = partStartIdx + 1;
    selectedParts = zeros(numberOfThreads,1);
    selectedParts(1) = partStartIdx;
    while selectedPartCount<numberOfThreads && partEndIdx <= numel(vocabDescriptions) 
        matchFlag = false;
        description = vocabDescriptions{partEndIdx};
        
        for partItr = 1:selectedPartCount
            description2 = vocabDescriptions{selectedParts(partItr)};
            matchingCost = InexactMatch(description, description2, maxDistance);
            if matchingCost <= threshold
                matchFlag = true;
                break;
            end
        end
        
        if ~matchFlag 
            selectedPartCount = selectedPartCount+1;
            selectedParts(selectedPartCount) = partEndIdx;
        end
        partEndIdx = partEndIdx+1;
    end
end

%> Name: InexactMatch
%>
%> Description: Given two coarse part descriptions, this function tries to
%> match them with the lowest cost possible. Please note that this is
%> essentially graph matching problem, which is NP-Complete. Using very high
%> dimension will result in a drastic performance degradation.
%> Cost of replacing a node: 1
%> Cost of re-positioning a node: distance/maxDistancePossible
%>
%> @param description Description of the first composition. It is of the
%> form: nodeId1 posX posY;
%>       nodeId2 posX posY; ...
%> @param description2 Description of the second composition.
%> @param maxDistance Maximum distance in the node positioning space.
%>
%> @retval lowestCost Minimum matching score.
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 06.05.2014
function [lowestCost] = InexactMatch(description, description2, maxDistance)
    global simMat;
    
    % Get both descriptions to the same size.
    firstDesSize = size(description,1);
    secDesSize = size(description2,1);
    if firstDesSize > secDesSize
        description2 = [description2; inf(firstDesSize-secDesSize,3)];
    elseif firstDesSize<secDesSize
        description = [description; inf(secDesSize-firstDesSize,3)];
    end
    numberOfChildren = max(firstDesSize, secDesSize);
    
    % Get row permutations of the first description.
    rows = perms(1:numberOfChildren);
    
    % Compare each permutation of rows of description to description2. The
    % one which gets the minimum cost is our match.
    lowestCost = inf;
    numberOfRows = size(rows,1);
    for permItr = 1:numberOfRows
        currentCost = 0;
        comparedDescription = description(rows(permItr,:),:);

        % Get valid rows to compare.
        validEdges = ~isinf(comparedDescription(:,1)) & ...
            ~isinf(description2(:,1));
        
        % Estimate node-node distances.
        for nodeItr = 1:numberOfChildren
            if validEdges(nodeItr)
                currentCost = currentCost + simMat(comparedDescription(nodeItr,1), ...
                                            description2(nodeItr,1));
            end
        end

        % Estimate edge-edge distances.
        currentCost = currentCost + sum(sqrt(sum((comparedDescription(validEdges,2:3) - ...
                     description2(validEdges,2:3)).^2,2)))/maxDistance;
        currentCost = currentCost + numberOfChildren - nnz(validEdges);
        
        % Assign lowest cost if current cost is smaller.
        if currentCost<lowestCost
            lowestCost = currentCost;
            if lowestCost == 0
                break;
            end
        end
    end
end
