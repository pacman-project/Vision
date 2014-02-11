%> Name: runSubdue
%>
%> Description: This function runs SUBDUE (self) with the input image and
%> parameters defined in options. 
%>
%> @param vocabLevel Input vocabulary level. If preDefinedSearch is 0,
%> ignored. If 1, compositions in this vocabulary level are detected in
%> graphLevel.
%> @param graphLevel The current object graphs' level. The graphLevel's
%> nodes should be sorted by their image id, and then rfId. 
%> @param options Program options.
%> @param preDefinedSearch If true, pre-defined subs in vocabLevel are found
%> in graphLevel. If false, graphLevel is given to unsupervised discovery, and
%> detected parts are returned in vocabLevel, while their realizations reside
%> in graphLevel.
%>
%> @retval nextVocabLevel Next vocabulary level ([] if pre-defined search
%> is on). 
%> @retval nextGraphLevel Next object graphs' level.
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 05.02.2014
function [nextVocabLevel, nextGraphLevel] = runSubdue(vocabLevel, graphLevel, options, preDefinedSearch)
    %% First thing we do is to convert vocabLevel and graphLevel into different data structures.
    % This process is done to assure fast, vectorized operations.
    % Initialize the priority queue.
    
    %% Some variables are global since they are accessed often in sub-functions.
    % Copying into function spaces and back is inefficient for this case.
    % They will be removed at the end of this function.
    global globalEdges globalRealNodeIds globalCenterSigns globalLabels globalRealNeighbors
    
    bestSubs = [];
    parentSubs = [];
    nextVocabLevel = [];
    nextGraphLevel = [];
    
    % Helper data structures.
    centerIdx = cat(1, graphLevel.isCenter) > 0;
    
    globalEdges = cat(1, graphLevel.adjInfo);
    
    % If no edges are present, time to return.
    if isempty(globalEdges)
        return;
    end
    
    globalEdges = globalEdges(:,1:3);
    globalRealNodeIds = cat(1, graphLevel.realNodeId);
    globalRealNeighbors = cat(1, {graphLevel(centerIdx).adjInfo});
%     for itr = 1:numel(globalRealNodeIds)
%         if size(globalRealNeighbors(itr),1) > 0
%             neighbors = globalRealNeighbors(itr);
%             globalRealNeighbors(itr) = globalRealNodeIds(neighbors(:,2));
%         else
%             globalRealNeighbors(itr) = [];
%         end
%     end
    nonemptyNeighborIdx = cellfun(@(x) ~isempty(x), globalRealNeighbors);
    globalRealNeighbors(nonemptyNeighborIdx) = cellfun(@(x) globalRealNodeIds(x(:,2)), ...
        globalRealNeighbors(nonemptyNeighborIdx), 'UniformOutput', false);
    globalCenterSigns = cat(1, graphLevel(centerIdx).sign);
    globalLabels = cat(1, graphLevel.labelId);
    
    %% Step 1:Find single-vertex subs, and put them into beamSubs.
    singleNodeSubs = getSingleNodeSubs(find(centerIdx), options);
    parentSubs = addToQueue(singleNodeSubs, parentSubs, options.subdue.beam);
    
    %% Step 2: Main loop
    startTime = tic;
    endFlag = 0;
    while ~isempty(parentSubs)
        extendedSubs = [];
        childSubArr = cell(numel(parentSubs),1);
        parentItr = 1;
        while ~isempty(parentSubs)
            %% Step 2.1: Get head of parentSubs.
            [parentSub, parentSubs] = getTopQueue(parentSubs);
            
            % Check time. If it took too long, end processing. 
            % Check parent's size. If it is too large, end processing.
            elapsedTime = toc(startTime);
            if elapsedTime > options.subdue.maxTime || ...
                    size(parentSub.edges,1) >= (options.subdue.maxSize-1)
                endFlag = 1;
                break;
            end
            %% Step 2.2: Extend head in all possible directions into childSubs.
            childSubArr(parentItr) = {extendSub(parentSub, options)};
            parentItr = parentItr + 1;
        end
        
        % If processing has been stopped due to time constraints, return bestSubs.
        if endFlag
           break; 
        end
        
        childSubs = cat(2, childSubArr{:});
        
        % If no children are generated, return.
        if numel(childSubs) < 1
           break; 
        end
        %% Step 2.3: Remove duplicates from childSubs.
        childSubs = removeDuplicateSubs(childSubs);
        
        %% Step 2.3: Evaluate childSubs, find their instances.
        childSubs = evaluateSubs(childSubs, options);

        %% Step 2.4: Add childSubs to extendedSubs.
        extendedSubs = addToQueue(childSubs, extendedSubs, options.subdue.beam);
        if numel(childSubs(1).edges) >= (options.subdue.minSize-1)
            bestSubs = addToQueue(childSubs, bestSubs, options.subdue.nsubs);    
        end
        
        %% Step 2.5: Swap expandedSubs with parentSubs.
        parentSubs = extendedSubs;
    end
    
    %% Step 3: Create nextVocabLevel and nextGraphLevel from bestSubs.
    numberOfBestSubs = numel(bestSubs);
    
    % If no subs detected, return.
    if numberOfBestSubs < 1
       vocabLevel = [];
       graphLevel = [];
       return;
    end
    
    %% Allocate space for new graphLevel and vocabLevel.
    instances = {bestSubs.instances};
    numberOfInstances = sum(cellfun(@(x) numel(x), instances));
    clear vocabLevel graphLevel
    vocabLevel(numberOfBestSubs) = options.vocabNode;
    graphLevel(numberOfInstances) = options.graphNode;
    
    %% Assign instances of each discovered substructure to graph node. 
    % We store substructures in vocabLevel array.
    for bestSubItr = 1:numberOfBestSubs
       % Assign label of sub.
       vocabLevel(bestSubItr).label = num2str(bestSubItr);
       
       % Assign sub edges (definition)
       numberOfEdges = size(bestSubs(bestSubItr).edges,1);
       if numberOfEdges > 0
           subEdges = [ones(numberOfEdges,1), ...
               (2:1:(numberOfEdges+1))', bestSubs(bestSubItr).edges(:,1), ones(numberOfEdges,1)];
           vocabLevel(bestSubItr).adjInfo = subEdges;
           vocabLevel(bestSubItr).children = [bestSubs(bestSubItr).centerId; bestSubs(bestSubItr).edges(:,2)]';
       else
           vocabLevel(bestSubItr).children = bestSubs(bestSubItr).centerId;
       end
    end
    
   %% Assign instances to graphLevel.
   % First, get labelId, children and children adjacency information
   % from the output.
   
   % Calculate labelIds array.
   labelIds = num2cell((1:numberOfBestSubs)');
   labelIds = cellfun(@(x,y) repmat(x, [numel(y),1]), labelIds, instances', 'UniformOutput', false);
   labelIds = cat(1, labelIds{:});
   labelIds = num2cell(labelIds);
   
   % Calculate children and the adjacency information linking children.
   edges = cellfun(@(x) {x.edges}, instances, 'UniformOutput', false);
   edges = [edges{:}];
   children = (cellfun(@(x) [x(1,1), x(x(:,5) > 0, 2)'], edges, 'UniformOutput', false))';
   childrenAdjInfo = cellfun(@(x) x(x(:,5) > 0, 1:3), edges, 'UniformOutput', false);
   childrenAdjInfo = cellfun(@(x) [x, ones([size(x,1), 1])], childrenAdjInfo', 'UniformOutput', false);
   
   assgnArr = [labelIds, children, childrenAdjInfo];

   % Out of this information, create graphLevel.
   [graphLevel.labelId, graphLevel.children, graphLevel.childrenAdjInfo] = deal(assgnArr{:});
       
   nextVocabLevel = vocabLevel;
   nextGraphLevel = graphLevel;
   %% TODO: Implement pre-defined search mode.
   
   
   
   
       
    clearvars -global
end

function singleNodeSubs = getSingleNodeSubs(centerIdx, options)
    global globalCenterSigns globalEdges globalLabels
    centerLabels = globalLabels(centerIdx,:);
    augmentedEdges = [globalEdges, globalLabels(globalEdges(:,2)), zeros(size(globalEdges,1), 1)];
    centerEdges = arrayfun(@(x) augmentedEdges(augmentedEdges(:,1) == x,:), centerIdx, 'UniformOutput', false);
    
    nodeAssgnArr = [num2cell([centerIdx, globalCenterSigns]), centerEdges];
    centerLabelIds = unique(centerLabels)';
    singleNodeSubs(numel(centerLabelIds)) = options.selfSubdue.sub;
    
    %% For each center node label type, we create a substructure.
    for centerLabelItr = 1:numel(centerLabelIds)
        selfCenterIdx = centerLabels == centerLabelIds(centerLabelItr);
        centerInstances = find(selfCenterIdx);
        
        %Assign center id.
        singleNodeSubs(centerLabelItr).centerId = centerLabelIds(centerLabelItr);
        
        % Give maximum score so that it is at the top of things.
        singleNodeSubs(centerLabelItr).mdlScore = numel(centerInstances);
        singleNodeInstances(numel(centerInstances)) = options.selfSubdue.instance;
        
        % Fill in instance information. 
        subNodeAssgnArr = nodeAssgnArr(selfCenterIdx,:);
        [singleNodeInstances.centerIdx, singleNodeInstances.sign, singleNodeInstances.edges] = deal(subNodeAssgnArr{:});
        
        singleNodeSubs(centerLabelItr).instances = singleNodeInstances;
        clear singleNodeInstances;
    end
end

function [extendedSubs] = extendSub(sub, options)
    % Get unused edges from sub's instances.
    allEdges = {sub.instances.edges};
    allUnusedEdges = cellfun(@(x) x(~x(:,5),:), allEdges, 'UniformOutput', false);
    unusedEdges = cat(1, allUnusedEdges{:});
    
    % Record which edge belongs to which instance. 
    allEdgeInstanceIds = zeros(size(unusedEdges,1),1);
    itrOffset = 1;
    unusedEdgeCount = cellfun(@(x) size(x,1), allUnusedEdges);
    for itr = 1:numel(allUnusedEdges)
        beginOffset = itrOffset;
        allEdgeInstanceIds(beginOffset:(beginOffset+unusedEdgeCount(itr))) = itr;
        itrOffset = itrOffset + unusedEdgeCount(itr);
    end
    
    % Get unique rows of [edgeLabel, secondVertexLabel]
    uniqueEdgeTypes = unique(unusedEdges(:, 3:4), 'rows');
    numberOfEdgeTypes = size(uniqueEdgeTypes,1);
    
    % Extend the definition in sub with each edge type in uniqueEdgeTypes.
    % In addition, we pick suitable instances, add this edge, and mark used
    % field of relevant instances.
    if numberOfEdgeTypes == 0
        extendedSubs = [];
        return;
    end
        
    extendedSubs(numberOfEdgeTypes) = options.selfSubdue.sub;
    for edgeTypeItr = 1:numberOfEdgeTypes
        % Assign sub-definition type and other info
        newSub = sub;
        newSub.mdlScore = 0;
        newSub.edges = [newSub.edges; uniqueEdgeTypes(edgeTypeItr,:)];
        
        %% Find instances of this new sub, and mark the new edges as 'used' in each of its subs.
        % This process is crucial for fast calculation of DL.
        edgesToExtend = ismember(unusedEdges(:,3:4), uniqueEdgeTypes(edgeTypeItr,:), 'rows');
        edgeInstances = unusedEdges(edgesToExtend,:);
        edgeInstanceIds = allEdgeInstanceIds(edgesToExtend,:);
        newInstances = sub.instances(edgeInstanceIds);
        
        % Each added edge actually means a new instance. 
        for instanceItr = 1:size(edgeInstances,1)
            addedEdgeIdx = ismember(newInstances(instanceItr).edges(:, 1:2), edgeInstances(instanceItr,1:2), 'rows');
            newInstances(instanceItr).edges(addedEdgeIdx, 5) = 1;
        end
        
        instanceEdges = {newInstances.edges}';
        validInstanceIdx = cellfun(@(x) ~isempty(x), instanceEdges);
        newInstances = newInstances(validInstanceIdx);
        %% Before assigning new instances to newSub, we eliminate duplicates.
        % Because of the instance generation process, an instance has many
        % parse trees if it has same type of edges between its nodes. So,
        % we may have the very same instance twice or many times in
        % newInstances list.
        if isempty(newInstances)
            instanceEdges = {newInstances.edges}';
            usedNodes = cell2mat(cellfun(@(x) x(x(:,5)>0, 2)', instanceEdges, 'UniformOutput', false));
            instanceIdentifier = [cat(1, newInstances.centerIdx), usedNodes];
            [~, uniqueInstanceIdx, ~] = unique(instanceIdentifier, 'rows');
            newInstances = newInstances(uniqueInstanceIdx);
        end
        
        newSub.instances = newInstances;
        extendedSubs(edgeTypeItr) = newSub;
    end
end

function [subs] = evaluateSubs(subs, options)
    global globalRealNeighbors globalRealNodeIds globalCenterSigns
    
    numberOfSubs = numel(subs);
    if strcmp(options.subdue.evalMetric, 'mdl')
        % This part is crucial! Essentially vbits + ebits + rbits defined for 
        % each neighbor node pair. The simplified graph structure allows for
        % such approximations of mdl score to work well. So, instead of
        % calculating DL on the fly, we calculate a score of combining
        % centerNode(n,1) with centerNode(n,2) by estimating compression amount
        % of removing centerNode(n,2) and adding its edges to centerNode(n,1).
        % Common neighbors need not be linked twice, and this is where
        % compression comes from. In addition, because the second node and the
        % edge linking first node to second node is removed, this brings
        % further compression in the main graph.
        % First term (2 * (pairwiseCommonNeighbors+1)): No need to hold common
        % neighbors' edges, in addition to the compressed node 
        % (each edge is represented with minimum 2 integers, edge label and neighbor id)
        % Second term (1): No need to hold compressed (neighbor vertex) any
        % more. Each vertex, essentially, consists of a label id, and
        % everything else is redundant information. 
        % Hint: Storage for different types of information is considered as
        % same for efficiency purposes (integer for all types, no compression
        % in label representation, or any other point).
        
        for subItr = 1:numberOfSubs
            instanceCenterSign = globalCenterSigns(globalRealNodeIds(cat(1, subs(subItr).instances.centerIdx)));
            instanceEdges = {subs(subItr).instances.edges};
            instanceSecondaryNeighbors = cellfun(@(x) globalRealNodeIds(x(x(:,5)>0,2)), instanceEdges, 'UniformOutput', false);
            instanceSecondaryNeighborIdx = cellfun(@(x) cat(1, globalRealNeighbors{x}), instanceSecondaryNeighbors, 'UniformOutput', false);
            instanceNeighborIdx = cellfun(@(x) globalRealNodeIds(x(:,2)), instanceEdges, 'UniformOutput', false);
            
            intersectingNodeCount = cellfun(@(x,y) numel(intersect(x,y)), instanceNeighborIdx, instanceSecondaryNeighborIdx);
            
            %% Estimate DL reduction.
            dlDiff = sum(intersectingNodeCount(instanceCenterSign == 1)) - sum(intersectingNodeCount(instanceCenterSign == 0));
            subs(subItr).mdlScore = dlDiff - size(subs(subItr).edges,1);
        end
    else
        for subItr = 1:numberOfSubs
            subs(subItr).mdlScore = size(subs(subItr).edges,1) * numel(subs(subItr).instances);
        end
    end
end

function childSubs = removeDuplicateSubs(childSubs)
    if isempty(childSubs)
       return; 
    end
    edges = {childSubs.edges};
    edges = cellfun(@(x) sortrows(x), edges, 'UniformOutput', false);
    edges = cellfun(@(x) (x(:)), edges, 'UniformOutput', false);
    edges = cell2mat(edges)';
    subIdentifiers = [cat(1, childSubs.centerId), edges];
    [~, uniqueSubIdentifiers, ~] = unique(subIdentifiers, 'rows', 'stable');
    childSubs = childSubs(uniqueSubIdentifiers);
end

%> Name: addToQueue
%>
%> Description: Adds subs in beamQueue, which is a priority queue that orders 
%> the substructures inside by their mdlScore.
%> 
%> @param subs An array of substructures to add to beamQueue.
%> @param queue The priority queue of substructures, ordered by
%> mdlScore. Best sub has highest score.
%> @param maxSize Maximum number of elements allowed in queue.
%>
%> @retval beamQueue The priority queue of substructures, ordered by
%> mdlScore.
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 06.02.2014
function [queue] = addToQueue(subs, queue, maxSize)
    addedQueue = [subs, queue];
    if maxSize>numel(addedQueue)
        maxSize = numel(addedQueue);
    end
    if isempty(addedQueue)
       queue = [];
       return; 
    end
    sortedQueue = nestedSortStruct(addedQueue, {'mdlScore'}, -1);
    queue = sortedQueue(1:maxSize);
end

%> Name: getTopQueue
%>
%> Description: Get top sub from queue, and delete it from queue.
%> 
%> @param queue Priority queue sorted by mdlScore.
%>
%> @retval topSub Best sub in restOfQueue.
%> @retval restOfQueue Rest of queue, minus topSub.
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 06.02.2014
function [topSub, restOfQueue] = getTopQueue(queue)
    if isempty(queue)
        restOfQueue = [];
        topSub = [];
    else
        topSub = queue(1);
        if numel(queue)>1
            restOfQueue = queue(2:end);
        else
            restOfQueue = [];
        end
    end
end


