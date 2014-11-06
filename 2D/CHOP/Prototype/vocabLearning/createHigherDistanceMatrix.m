function [vocabLevel, graphLevel, newDistanceMatrix, subClasses] = createHigherDistanceMatrix(vocabLevel, nodeDistanceMatrix, prevGraphSize, prevGraphData, options)
    setSize = 20; % Parallelization parameter. Probably sole hard coded parameter in the system. It is a thing of beauty.
    regularizationParam = (options.subdue.maxSize * 2) - 1; % Maximum size of a part (n nodes + n-1 edges)
    threshold = options.subdue.threshold * regularizationParam; % Hard threshold for cost of matching two subs.
    subClasses = [];
    newDistanceMatrix = [];
    
    % Get program params.
    evalMetric = options.subdue.evalMetric;
    mdlNodeWeight = options.subdue.mdlNodeWeight;
    mdlEdgeWeight = options.subdue.mdlEdgeWeight;
    isMDLExact = options.subdue.isMDLExact;
    overlap = options.subdue.overlap;
    isSupervised = options.subdue.supervised;
    
    %% If no current modes exist, we eliminate them based on their children only.
    if isempty(currentModes)
        vocabDescriptions = {vocabLevel.children};
        vocabDescriptions = cellfun(@(x) mat2str(sort(x)), vocabDescriptions, 'UniformOutput', false);
        [~, IA, IC] = unique(vocabDescriptions, 'stable');
        clear vocabDescriptions;
        vocabLevel = vocabLevel(1, IA);

        % Update graph level.
        labelIds = cat(1, graphLevel.labelId);
        IC = num2cell(int32(IC(labelIds)));
        [graphLevel.labelId] = deal(IC{:});
        clear IC;
    else
        %% If current modes do exist, we can use the positioning information to eliminate parts further.
        currentModes = currentModes(:,3:4);
        currentModes((size(currentModes,1)+1),:) = [0, 0];

        % Get unique sub descriptions independent of the center node. Each
        % description is a string that encodes coarse structure of the sub.
        vocabNodeLabels = {vocabLevel.children};
        vocabNodeLabels = cellfun(@(x) double(x), vocabNodeLabels, 'UniformOutput', false);
        vocabEdges = {vocabLevel.adjInfo};
        vocabEdges = cellfun(@(x) double(x), vocabEdges, 'UniformOutput', false);
        newMode = size(currentModes,1);
        vocabNeighborModes = cellfun(@(x) [newMode; x(:,3)], vocabEdges, 'UniformOutput', false);
        vocabNodePositions = cellfun(@(x) currentModes(x,:) - repmat(min(currentModes(x,:)), numel(x), 1), vocabNeighborModes, 'UniformOutput', false);

        % Maximum distance between two nodes in 2D space.
        maxDistance = max(sqrt(sum(currentModes.^2,2)))*2;
        vocabSortOrder = cell(size(vocabLevel,1),1);
        for vocabNodeItr = 1:numel(vocabLevel)
            [~, vocabSortOrder{vocabNodeItr}] = sortrows(vocabNodePositions{vocabNodeItr});
        end
        vocabDescriptions = cellfun(@(x,y,z) [x(z)', y(z,:)], vocabNodeLabels, vocabNodePositions, vocabSortOrder, 'UniformOutput', false);
        clear vocabNodeLabels vocabEdge svocabNeighborModes vocabNodePositions vocabSortOrder;

        %% For each sub, we will find matching ones with a cost.
        numberOfSubs = numel(vocabDescriptions);
        matchedSubs = zeros(1, numberOfSubs)>0;
        subClasses = int32((1:numberOfSubs)');
        newDistanceMatrix = zeros(numberOfSubs);

        % Here, we parallelize the part comparison process. In order to
        % do that, numberOfThreads best parts are selected, while the
        % rest of the parts are compared with the selected ones. When
        % comparison is over, the data structures are updated, and a
        % new set of parts are selected from remaining list.
        partStartIdx = 1;
        oldDistMat = distanceMatrix;

        while partStartIdx < numberOfSubs
            % Select first numberOfThreads parts.
            [selectedParts, partEndIdx, ~, ~, ~] = GetBestParts(vocabDescriptions, partStartIdx, ...
                matchedSubs, subClasses, setSize, maxDistance, distanceMatrix, threshold);
%             if numel(selectedParts)>1
%                 bestEntries = [bestEntries, zeros(size(bestEntries,1),1)]; %#ok<AGROW>
%                 bestEntries = [zeros(1, size(bestEntries,2)); bestEntries]; %#ok<AGROW>
%                 bestEntries = bestEntries + bestEntries';
%                 newDistanceMatrix(selectedParts, selectedParts) = bestEntries;
%             end
            
%             if partEndIdx>numberOfSubs
%                 break;
%             end

            distanceMatrixEntries = cell(numberOfSubs,1);
            % If set is empty, we're at the end of the list.
            if isempty(selectedParts)
                break;
            end

            % Check against the rest of the parts.
            parfor partItr = (partStartIdx+1):numberOfSubs
                if ~matchedSubs(partItr)
                    description = vocabDescriptions{partItr};
                    newEntries = zeros(1, numel(selectedParts));
                    for partItr2 = 1:numel(selectedParts)
                        description2 = vocabDescriptions{selectedParts(partItr2)}; %#ok<PFBNS>
                        matchingCost = InexactMatch(description, description2, maxDistance, oldDistMat, threshold);
                        if matchingCost <= threshold && selectedParts(partItr2) < partItr
                            matchedSubs(partItr) = 1;
                            subClasses(partItr) = selectedParts(partItr2);
                            break; 
                        else
%                             if matchingCost~=newDistanceMatrix(selectedParts(partItr2), partItr)
%                                 1
%                             end
                            newEntries(partItr2) = matchingCost;
                        end
                    end
                    if ~matchedSubs(partItr)
                        distanceMatrixEntries{partItr} = newEntries;
                    end
                end
            end

            % Update the similarity matrix based on the new entries.
            allEntries = cat(1, distanceMatrixEntries{:});
            nonemptyIdx = cellfun(@(x) ~isempty(x), distanceMatrixEntries);
            newDistanceMatrix(nonemptyIdx, selectedParts) = allEntries;
            newDistanceMatrix(selectedParts, nonemptyIdx) = allEntries';
            partStartIdx = partEndIdx;
            clear distanceMatrixEntries;
        end

        % Process the last sub, if it has not been assigned to any
        % other sub.
        if ~matchedSubs(end)
            subClasses(end) = numberOfSubs;
        end

        % Create new distance matrix.
        newDistanceMatrix = newDistanceMatrix/max(max(newDistanceMatrix));
        
        %% We calculate scores for each vocabulary node and re-order them again after combining parts.
        % Calculate scores for each node.
        allEdges = prevGraphData.allEdges;
        allEdgeNodePairs = prevGraphData.allEdgeNodePairs;
        allSigns = prevGraphData.allSigns;
        bestSubs = prevGraphData.bestSubs;
        bestSubs = bestSubs(1:numel(matchedSubs));        
        % Make the instance lists uniform so that they can be concatenated
        % easily.
        for subItr = 1:numel(bestSubs)
            if size(bestSubs(subItr).instances,1) ~= 1
                bestSubs(subItr).instances = bestSubs(subItr).instances';
            end
        end
        
        for subItr = 1:numel(matchedSubs)
            if ~matchedSubs(subItr)
                bestSubs(subItr).instances = cat(2, bestSubs(subClasses == subItr).instances);
                if isSupervised
                    categoryArr = double([bestSubs(subItr).instances.category]);
                    weight = nnz(categoryArr == mode(categoryArr)) / numel(categoryArr);
                else
                    weight = 1;
                end
                newMDLScore = weight * getSubScore(bestSubs(subItr), allEdges, allEdgeNodePairs, evalMetric, ...
                allSigns, mdlNodeWeight, mdlEdgeWeight, overlap, isMDLExact);
                vocabLevel(subItr).mdlScore = newMDLScore;
                if strcmp(evalMetric, 'mdl')
                    vocabLevel(subItr).normMdlScore = 1 - (newMDLScore / prevGraphSize);
                end
            else
                vocabLevel(subItr).mdlScore = -inf;
                vocabLevel(subItr).normMdlScore = inf;
            end
        end
        
        % Sort the new parts based on their mdl scores.
        mdlScores = [vocabLevel.mdlScore];
        [~, sortedIdx] = sort(mdlScores, 'descend');
        vocabLevel = vocabLevel(sortedIdx);
        sortedAssgnArr = int32(1:numel(vocabLevel));
        sortedAssgnArr(sortedIdx) = sortedAssgnArr;
        
        % Assign labels to each vocabulary node.
        for subItr = 1:numel(vocabLevel)
            vocabLevel(subItr).label = int32(subItr);
        end
        
        % Link each node in the redundant vocabulary to the updated nodes, 
        % too.
        subClasses = subClasses(sortedIdx);
        subClasses = sortedAssgnArr(subClasses)';
        
        % Update remaining data structures for output.
        newDistanceMatrix = newDistanceMatrix(sortedIdx, sortedIdx);
        graphLabelIds = [graphLevel.labelId];
        updatedGraphLabelIds = sortedAssgnArr(graphLabelIds);
        updatedGraphLabelIds = num2cell(updatedGraphLabelIds);
        [graphLevel.labelId] = deal(updatedGraphLabelIds{:});
        clear updatedGraphLabelIds;
        
        %% Eliminate parts which have negative mdl values. This design
        % choice helps us to get rid of parts which actually do not
        % compress the object graphs, but rather increase its size.
        mdlScores = [vocabLevel.mdlScore];
        validMDLScoreIdx = mdlScores >= 0 | mdlScores == -inf;
        invalidMDLScoreIdx = subClasses(~validMDLScoreIdx);
        validMDLScoreIdx = ~ismember(subClasses, invalidMDLScoreIdx);
        % Update relevant data structures.
        vocabLevel = vocabLevel(validMDLScoreIdx);
        % Assign labels to each vocabulary node.
        for subItr = 1:numel(vocabLevel)
            vocabLevel(subItr).label = int32(subItr);
        end
        subClasses = subClasses(validMDLScoreIdx);
        newDistanceMatrix = newDistanceMatrix(validMDLScoreIdx,validMDLScoreIdx);
        
        % Update graph label ids so that the remaining ids are consecutive 
        % and they do not contain any gaps (deleted parts). Redundant part
        % ids are also removed from the object graph (graphLevel), replaced
        % with real labels.
        if isempty(graphLevel)
           return; 
        end
        graphLabelIds = [graphLevel.labelId];
        remainingComps = find(validMDLScoreIdx);
        validGraphLevelIdx = ismember(graphLabelIds, remainingComps);
        graphLevel = graphLevel(validGraphLevelIdx);
        graphLabelIds = graphLabelIds(validGraphLevelIdx);
        sortedAssgnArr = zeros(max(remainingComps),1);
        sortedAssgnArr(remainingComps) = 1:numel(remainingComps);
        graphLabelIds = subClasses(sortedAssgnArr(graphLabelIds));
        graphLabelIdsCell = num2cell(graphLabelIds);
        if ~isempty(graphLevel)
            [graphLevel.labelId] = deal(graphLabelIdsCell{:});
        else
            return;
        end
        
        %% Update graph level's ids, and sort them based on labelIds.
        graphLabelIds = cat(1, graphLevel.labelId);
        imageIds = cat(1, graphLevel.imageId);
        arrayToSort = [imageIds, graphLabelIds];
        [~, sortedIdx2] = sortrows(arrayToSort);
        graphLevel = graphLevel(sortedIdx2);
    end    
end