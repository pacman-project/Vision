%> Name: runSubdue
%>
%> Description: This function runs SUBDUE (self) with the input image and
%> parameters defined in options. 
%>
%> @param vocabLevel Input vocabulary level. If preDefinedSearch is 0,
%> ignored. If 1, compositions in this vocabulary level are detected in
%> graphLevel.
%> @param graphLevel The current object graphs' level. The graphLevel's
%> nodes should be sorted by their image id. 
%> @param oppositeModes List of the form [oppMode1; oppMode2; ...].
%> label1 encodes exact opposite geometric information of each mode in its 
%> index, and same applies to every row of this list. If empty, simply ignored.
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
function [nextVocabLevel, nextGraphLevel] = runSubdue(vocabLevel, graphLevel, oppositeModes, options, ~)
    %% First thing we do is to convert vocabLevel and graphLevel into different data structures.
    % This process is done to assure fast, vectorized operations.
    % Initialize the priority queue.
    
    % Some variables are global since they are accessed often in sub-functions.
    % Copying into function spaces and back is inefficient for this case.
    % They will be removed at the end of this function.
    global globalEdges globalLabels globalNeighbors globalSigns
    
    bestSubs = [];
    parentSubs = [];
    nextVocabLevel = [];
    nextGraphLevel = [];
    
    % Helper data structures.
    globalEdges = {graphLevel.adjInfo}';
    
    % If no edges are present, time to return.
    if isempty(globalEdges)
        return;
    end
    globalNeighbors = globalEdges;
    nonemptyEdgeIdx = cellfun(@(x) ~isempty(x), globalEdges);
    globalNeighbors(nonemptyEdgeIdx) = cellfun(@(x) x(:,2), ...
        globalEdges(nonemptyEdgeIdx), 'UniformOutput', false);
    globalLabels = cat(1, graphLevel.labelId);
    globalEdges(nonemptyEdgeIdx) = cellfun(@(x) [x(:,1:3), globalLabels(x(:,2))], ...
        globalEdges(nonemptyEdgeIdx), 'UniformOutput', false);
    globalSigns = cat(1, graphLevel.sign);
    
    %% Step 1:Find single-vertex subs, and put them into beamSubs.
    singleNodeSubs = getSingleNodeSubs(options);
    parentSubs = addToQueue(singleNodeSubs, parentSubs, options.subdue.beam);
    
    %% Step 2: Main loop
    startTime = tic;
    endFlag = 0;
    while ~isempty(parentSubs)
        extendedSubs = [];
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
            childSubs = extendSub(parentSub, options);
            
            % If no children are generated or time is up, return.
            if endFlag || numel(childSubs) < 1
               break; 
            end
            
            %% Step 2.3: Remove duplicates from childSubs.
            % Here, we check if childSubs has duplicates in extendedSubs.
            childSubs = removeDuplicateSubs([extendedSubs, childSubs]);
            childSubs = childSubs([childSubs.mdlScore] == 0);
            
            %% Step 2.4: Evaluate childSubs, find their instances.
            childSubs = evaluateSubs(childSubs, options);
            
            if isempty(childSubs)
                continue;
            end
            %% Step 2.5: Add childSubs to extendedSubs and bestSubs.
            extendedSubs = addToQueue(childSubs, extendedSubs, options.subdue.beam);
            if numel(childSubs(1).edges) >= (options.subdue.minSize-1)
                % Here, a check ensures that the subs to put in bestSubs
                % are distinct. Simply, an instance of any sub in bestSubs
                % should not have the exact same node list (children) as
                % any other sub in the list.
                [~,sortedIdx]=sort([childSubs.mdlScore]);
                childSubs=childSubs(sortedIdx);
                childSubs = removeEncodedSubs(childSubs, bestSubs, oppositeModes);
                
                % Add remaining child subs to best subs.
                bestSubs = addToQueue(childSubs, bestSubs, options.subdue.nsubs);    
            end
        end
            
        %% Step 2.6: Swap expandedSubs with parentSubs.
        parentSubs = extendedSubs;
    end
    
    %% Step 3: Create nextVocabLevel and nextGraphLevel from bestSubs.
    numberOfBestSubs = numel(bestSubs);
    
    % If no subs detected, return.
    if numberOfBestSubs < 1
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
    instanceOffset = 1;
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
       %% Assign instances to this sub.
       instances = bestSubs(bestSubItr).instances;
       numberOfInstances = numel(bestSubs(bestSubItr).instances);
       labelIds = num2cell(repmat(bestSubItr, numberOfInstances,1));
       
       %% Create required fields such as centerIdx, edges, children and sign array to assign to instances.
       centerIdx = [instances.centerIdx];
       centerIdxCellArr = num2cell(centerIdx');
       numberOfInstances = numel(centerIdx);
       instanceEndOffset = instanceOffset + numberOfInstances - 1;
       edges = globalEdges(centerIdx);
       edgeIdx = {instances.edges}';
       instanceEdges = cellfun(@(x,y) x(y,:), edges, edgeIdx, 'UniformOutput', false);
       instanceChildren = cellfun(@(x,y) [x, y(:,2)'], centerIdxCellArr, instanceEdges, 'UniformOutput', false);
       instanceSigns = num2cell(globalSigns(centerIdx));
       
       [graphLevel(instanceOffset:instanceEndOffset).labelId] = deal(labelIds{:});
       [graphLevel(instanceOffset:instanceEndOffset).children] = deal(instanceChildren{:});
       [graphLevel(instanceOffset:instanceEndOffset).sign] = deal(instanceSigns{:});
       instanceOffset = instanceOffset + numberOfInstances;
       
       
       clear centerIdx centerIdxCellArr edges instanceEdges instanceChildren instanceSigns
    end
    
   nextVocabLevel = vocabLevel;
   nextGraphLevel = graphLevel;
   %% TODO: Implement pre-defined search mode.
   
   
    clearvars -global
end


%> Name: getSingleNodeSubs
%>
%> Description: getSingleNodeSubs is used to obtain single-node subs from
%> labels of all instances in globalLabels. The result is a number of
%> subs representing compositions each having their instances.
%> 
%> @param options Program options.
%>
%> @retval singleNodeSubs Substructure list of single-node subs.
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 24.02.2014
function singleNodeSubs = getSingleNodeSubs(~)
    global globalLabels globalSigns
    numberOfSubs = max(globalLabels);
    singleNodeSubs(numberOfSubs) = Substructure();
    
    %% For each center node label type, we create a substructure.
    for subItr = 1:numberOfSubs
        subCenterIdx = globalLabels == subItr;
        instances = find(subCenterIdx);
        numberOfInstances = numel(instances);
        
        %Assign center id.
        singleNodeSubs(subItr).centerId = subItr;
        
        % Give maximum score so that it is at the top of things.
        singleNodeSubs(subItr).mdlScore = numberOfInstances;
        singleNodeInstances(numberOfInstances) = Instance();
        
        % Fill in instance information. 
        instanceIdx = globalLabels == subItr;
        subNodeAssgnArr = num2cell([find(instanceIdx), globalSigns(instanceIdx,1)]);
        [singleNodeInstances.centerIdx, singleNodeInstances.sign] = deal(subNodeAssgnArr{:});
        
        singleNodeSubs(subItr).instances = singleNodeInstances;
        clear singleNodeInstances;
    end
end

%> Name: extendSub
%>
%> Description: extendSub(..) extends 'sub' in all possible ways
%> by extending its instances, and returns the new sub list 'extendedSubs'
%> along with instances of each sub in returned list.
%> 
%> @param sub Sub to be extended.
%> @param options Program options.
%>
%> @retval extendedSubs Extended sub list.
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 24.02.2014
function [extendedSubs] = extendSub(sub, options)
    global globalEdges
    
    centerIdx = cat(1,sub.instances.centerIdx);
    subAllEdges = globalEdges(centerIdx);
    % Get unused edges from sub's instances.
    allUsedEdgeIdx = {sub.instances.edges}';
    allUnusedEdges = cellfun(@(x,y) x(setdiff(1:size(x,1),y),:), subAllEdges, allUsedEdgeIdx, 'UniformOutput', false);
    allUnusedInstanceEdges = allUnusedEdges;
    
    % Record which edge belongs to which instance. 
    allEdgeInstanceIds = zeros(size(allUnusedEdges,1),1);
    itrOffset = 1;
    unusedEdgeCount = cellfun(@(x) size(x,1), allUnusedEdges);
    for itr = 1:numel(allUnusedEdges)
        beginOffset = itrOffset;
        allEdgeInstanceIds(beginOffset:(beginOffset+unusedEdgeCount(itr))) = itr;
        itrOffset = itrOffset + unusedEdgeCount(itr);
    end         
    allUnusedEdges = cat(1, allUnusedEdges{:});
    
    % Get unique rows of [edgeLabel, secondVertexLabel]
    uniqueEdgeTypes = unique(allUnusedEdges(:, 3:4), 'rows');
    
    % Discard any edge types already existing in sub definition.
    if ~isempty(sub.edges)
        uniqueEdgeTypes = setdiff(uniqueEdgeTypes, sub.edges, 'rows');
    end
    
    % Extend the definition in sub with each edge type in uniqueEdgeTypes.
    % In addition, we pick suitable instances, add this edge, and mark used
    % field of relevant instances.
    numberOfEdgeTypes = size(uniqueEdgeTypes,1);
    if numberOfEdgeTypes == 0
        extendedSubs = [];
        return;
    end
        
    extendedSubs(numberOfEdgeTypes) = Substructure;
    for edgeTypeItr = 1:numberOfEdgeTypes
        % Assign sub-definition type and other info
        newSub = sub;
        newSub.mdlScore = 0;
        newSub.edges = [newSub.edges; uniqueEdgeTypes(edgeTypeItr,:)];
        
        %% Find instances of this new sub, and mark the new edges as 'used' in each of its subs.
        % This process is crucial for fast calculation of DL.
        edgesToExtendIdx = allUnusedEdges(:,3) == uniqueEdgeTypes(edgeTypeItr,1) & ...
            allUnusedEdges(:,4) == uniqueEdgeTypes(edgeTypeItr,2);
        edgesToExtend = allUnusedEdges(edgesToExtendIdx,:);
        edgeInstanceIds = allEdgeInstanceIds(edgesToExtendIdx,:)';
        newInstances = sub.instances(edgeInstanceIds);
        
        % Each added edge actually means a new instance. 
        for instanceItr = 1:numel(edgeInstanceIds)
            instanceEdges = allUnusedInstanceEdges{edgeInstanceIds(instanceItr)};
            if isempty(instanceEdges)
               continue; 
            end
            addedEdgeIdx = find(instanceEdges(:, 2) == edgesToExtend(instanceItr,2));
            newInstances(instanceItr).edges = [newInstances(instanceItr).edges; addedEdgeIdx];
        end
        
%        % Eliminate those instances not having any edges.
%         instanceEdges = {newInstances.edges}';
%         validInstanceIdx = cellfun(@(x) ~isempty(x), instanceEdges);
%         newInstances = newInstances(validInstanceIdx);
        %% Before assigning new instances to newSub, we eliminate duplicates.
        % Because of the instance generation process, an instance has many
        % parse trees if it has same type of edges between its nodes. So,
        % we may have the very same instance twice or many times in
        % newInstances list.
%         if ~isempty(newInstances)
%             allNewInstanceEdges = allUnusedInstanceEdges(edgeInstanceIds);
%             instanceEdges = {newInstances.edges}';
%             usedNodes = cell2mat(cellfun(@(x,y) x(y, 2)', allNewInstanceEdges, instanceEdges, 'UniformOutput', false));
%             instanceIdentifier = [cat(1, newInstances.centerIdx), usedNodes];
%             [~, uniqueInstanceIdx, ~] = unique(instanceIdentifier, 'rows');
%             newInstances = newInstances(uniqueInstanceIdx);
%         end
        
        newSub.instances = newInstances;
        extendedSubs(edgeTypeItr) = newSub;
    end
end

%> Name: evaluateSubs
%>
%> Description: The evaluation function of SUBDUE. Based on the
%> evaluation metric, the value of each substructure in subs list is
%> calculated, and saved within subs. The MDL calculation takes place here.
%> 
%> @param subs Sub list which will be evaluated.
%> @param options Program options.
%>
%> @retval subs Evaluated sub list.
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 24.02.2014
function [subs] = evaluateSubs(subs, options)
    global globalSigns globalEdges
    
    numberOfSubs = numel(subs);
    if strcmp(options.subdue.evalMetric, 'mdl')
        
        % This part is crucial! Essentially vbits + ebits + rbits defined for 
        % each neighbor node pair. The simplified graph structure allows for
        % easy calculation of mdl score. We calculate a score of combining
        % centerNode(n,1) with centerNode(n,2) by estimating compression amount
        % of removing centerNode(n,2) and adding its edges to centerNode(n,1).
        % Common neighbors need not be linked twice, and this is where
        % compression comes from. In addition, because the second node and the
        % edge linking first node to second node is removed, this brings
        % further compression in the main graph.
        % Hint: Storage for different types of information is considered as
        % same for efficiency purposes (integer for all types, no compression
        % in label representation, or any other point).
        for subItr = 1:numberOfSubs
            % Read signs and edges of the instances.
            centerIdx = cat(1, subs(subItr).instances.centerIdx);
            instanceEdges = globalEdges(centerIdx);
            instanceSigns = globalSigns(cat(1, subs(subItr).instances.centerIdx));
            instanceUsedEdgeIdx = {subs(subItr).instances.edges}';
            
            % Calculate direct neighbors and neighbors of neighbors for mdl
            % calculation.
            instanceAllNeighbors = cellfun(@(x) x(:,2), instanceEdges, 'UniformOutput', false);
            instanceNeighbors = cellfun(@(x,y) x(y,2), instanceEdges, instanceUsedEdgeIdx, 'UniformOutput', false);
            instanceNeighborEdges = cellfun(@(x) cat(1, globalEdges{x}), instanceNeighbors, 'UniformOutput', false);
            nonemptyNeighborIdx = ~cellfun(@(x) isempty(x), instanceNeighborEdges);
            instanceSecondaryNeighbors = cell(size(nonemptyNeighborIdx,1),1);
            instanceSecondaryNeighbors(nonemptyNeighborIdx) = cellfun(@(x) x(:,2), instanceNeighborEdges(nonemptyNeighborIdx), 'UniformOutput', false);
            
            % Calculate DL reduction for each instance. 
            % (numel(fastintersect(x,y)) : Number of edges to be removed to
            % represent secondary nodes' edges with the center node.
            % (numel(y) - numel(unique(y))) : Number of duplicates within
            % y. Each secondary neighbor (neighbor of neighbor) needs to be
            % represented only once. Constant terms have been discarded
            % from this equation.
            dlReductions = cellfun(@(x,y) numel(fastintersect(x,y)) + (numel(y) - numel(fast_unique(y))), instanceAllNeighbors, instanceSecondaryNeighbors);
            
            %% Estimate DL reduction.
            dlDiff = sum(dlReductions(instanceSigns == 1)) - sum(dlReductions(instanceSigns == 0));
            subs(subItr).mdlScore = dlDiff - size(subs(subItr).edges,1);
        end
    else
        for subItr = 1:numberOfSubs
            subs(subItr).mdlScore = size(subs(subItr).edges,1) * (numel(subs(subItr).instances)-1);
        end
    end
end

%> Name: removeDuplicateSubs
%>
%> Description: Given a substructure list, this function removes duplicates
%> from it by calculating a unique description for each sub, and then
%> checking for redundant descriptions (a.k.a. duplicate subs).
%> 
%> @param childSubs The substructure list from which duplicates will be
%> removed.
%>
%> @retval childSubs Unique substructures.
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 24.02.2014
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

%> Name: removeEncodedSubs
%>
%> Description: Given subs to check (subs2Check) and already considered
%> subs (encodedSubs), this function removes the subs from subs2Check and
%> returns remaining subs in validSubs.
%> 
%> @param subs2Check Subs each of which is considered for elimination if
%> they already exist either in encodedSubs, or subs2Check as duplicates.
%> @param encodedSubs Subs already encoded.
%> @param oppositeEdgeLabelList nx1 array indicating an edge type's reverse
%> geometric information if it is encoded by another edge. 
%> oppositeEdgeLabels(x)=y => oppositeEdgeLabels(y)=x.
%>
%> @retval validSubs Unique subs from subs2Check.
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 15.02.2014
function validSubs = removeEncodedSubs(subs2Check, encodedSubs, oppositeEdgeLabelList)
    validSubs = [];
    if isempty(subs2Check) 
        return;
    end
    %% First, eliminate duplicate subs in subs2Check.
    % Check number of nodes in subs2Check. It is assumed that all subs in
    % subs2Check have the same size. If size ~= 2, there is no ambiguity in
    % the geometry, and each sub encodes a different structure. No need to
    % check.
    numberOfNodes = size(subs2Check(1).edges,1) + 1;
    
    if numberOfNodes == 2
        % First, we check for duplicate subs in subs2Check.
        centerIds = cat(1, subs2Check.centerId);
        edges = {subs2Check.edges};
        descriptions = [centerIds, cat(1, edges{:})];
        reverseDescriptions = [descriptions(:,3), oppositeEdgeLabelList(descriptions(:,2)), descriptions(:,1)];

        % Get description matches.
        [~, reverseIdx] = ismember(reverseDescriptions, descriptions, 'rows');

        % Eliminate half of the matches, since one of the two matching subs
        % will be left in the sub list.
        reverseCheckArr = (1:numel(reverseIdx))';
        remainingSubIdx = reverseIdx <= reverseCheckArr;
        subs2Check = subs2Check(remainingSubIdx);
    else
        remainingSubIdx = ones(numel(subs2Check),1) > 0;
    end
    validSubs = subs2Check;
    
    if isempty(encodedSubs)
        return;
    end
    %% Second task involves checking subs2Check against encodedSubs.
    encodedSubsNodeCountArr = {encodedSubs.edges};
    encodedSubsNodeCountArr = cellfun(@(x) size(x,1), encodedSubsNodeCountArr)' + 1;
    encodedSubs = encodedSubs(encodedSubsNodeCountArr == numberOfNodes);
    
    if numberOfNodes == 2
        % Eliminated those not having exact same number of nodes. If empty,
        % return.
        if isempty(encodedSubs)
            return
        end
        
        % Get reverse descriptions of each sub in subs2Check.
        reverseDescriptions = reverseDescriptions(remainingSubIdx, :);

        % Get descriptions of each sub in encodedSubs.
        encodedCenterIds = cat(1, encodedSubs.centerId);
        encodedEdges = {encodedSubs.edges};
        encodedDescriptions = [encodedCenterIds, cat(1, encodedEdges{:})];

        % Eliminate those in subs2Check that match already encoded structure in
        % encodedSubs.
        remainingSubIdx = ~ismember(reverseDescriptions, encodedDescriptions, 'rows');
        subs2Check = subs2Check(remainingSubIdx);
    end
    validSubs = subs2Check;
    
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
    [~,sortedIdx]=sort([addedQueue.mdlScore], 'descend');
    sortedQueue=addedQueue(sortedIdx);
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

