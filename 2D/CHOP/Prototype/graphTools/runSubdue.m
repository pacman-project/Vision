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
%> label1 encodes exact opposite geometric information of each mode in its 
%> index, and same applies to every row of this list. If empty, simply ignored.
%> @param options Program options.
%>
%> Explanation of the algorithm:
%> 
%> Given Graph V = (N, E), with node list N and edge Eist E, 
%>       Queue = [], ExtendQueue = [], FinalQueue = [], options.
%> Queue, ExtendQueue, FinalQueue are ordered by mdl scores. They are
%> automatically trimmed to options.beam elements after any addition.
%> 1- Find single node subs in N and add them to Queue. Occurences of each
%> node type is counted, and Queue is ordered by the number of occurences.
%> 
%> 2- For each sub S in Queue,
%>      2.1- Extend S by exploring all possible extensions in G 
%>      into ChildSubs.
%> 
%>      2.2- Remove any duplicate subs from ChildSubs. 
%> 
%>      2.3- Evauate each sub in ChildSubs and assign their mdl scores.
%>
%>      2.4- Add each sub in ChildSubs into ExtendQueue.
%>
%> 3- Add each sub in ExtendQueue to FinalQueue. No duplicate subs are
%> added to FinalQueue if they already exist.
%> 
%> 4- Write final subs and their instances to vocabLevel, graphLevel.
%>
%> @retval nextVocabLevel Next vocabulary level ([] if pre-defined search
%> is on). 
%> @retval nextGraphLevel Next object graphs' level.
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 05.02.2014
%> Ver 1.1 on 01.09.2014 Removal of global parameters.
function [nextVocabLevel, nextGraphLevel] = runSubdue(vocabLevel, graphLevel, options)
    %% First thing we do is to convert vocabLevel and graphLevel into different data structures.
    % This process is done to assure fast, vectorized operations.
    % Initialize the priority queue.
    
    bestSubs = [];
    parentSubs = [];
    nextVocabLevel = [];
    nextGraphLevel = [];
    
    %% Get the parameters.
    evalMetric = options.subdue.evalMetric;
    beam = options.subdue.beam;
    nsubs = options.subdue.nsubs;
    maxTime = options.subdue.maxTime;
    maxSize = options.subdue.maxSize;
    numberOfThreads = options.numberOfThreads;
    
    %% Initialize data structures.
    display('[SUBDUE] Initializing data structures for internal use..');
    % Helper data structures.
    allEdges = {graphLevel.adjInfo}';
    
    % If no edges are present, time to return.
    if isempty(allEdges)
        return;
    end
    nonemptyEdgeIdx = cellfun(@(x) ~isempty(x), allEdges);
    allLabels = cat(1, graphLevel.labelId);
    allEdges(nonemptyEdgeIdx) = cellfun(@(x) [x(:,1:3), allLabels(x(:,2))], ...
        allEdges(nonemptyEdgeIdx), 'UniformOutput', false);
    allSigns = cat(1, graphLevel.sign);
    
    % Graph size formulation is very simple: 2 * #edges + #nodes. edges are
    % linked to their node1, so only need to keep node2 and edgeLabel
    % (2xint). Nodes only keep nodeLabel (1xint). We assume nearly all
    % edges are doubly linked, therefore leading to an approximate yet fast
    % calculation of graph size.
    graphSize = sum(cellfun(@(x) size(x,1), allEdges)) * 2 + numel(graphLevel);
    
    %% Step 1:Find single-vertex subs, and put them into beamSubs.
    display('[SUBDUE] Creating single node subs..');
    singleNodeSubs = getSingleNodeSubs(allLabels, allSigns);
    parentSubs = addToQueue(singleNodeSubs, parentSubs, options.subdue.beam);
    
    %% Step 2: Main loop
    startTime = tic;
    
    while ~isempty(parentSubs)
        % Check time. If it took too long, end processing. 
        % Check parent's size. If it is too large, end processing.
        elapsedTime = toc(startTime);
        if size(parentSubs(1).edges,1) >= (maxSize-1)
            display(['[SUBDUE] Maximum size of ' num2str(maxSize) ' nodes reached, terminating search.']);
            break;
        end
        if elapsedTime > maxTime 
            display(['[SUBDUE] Time is up! Breaking from process just before extending subs of size ' num2str(size(parentSubs(1).edges,1)+1) '..']);
            break;
        end
        
        %% All good, execution.
        display(['[SUBDUE/Parallel] Expanding subs of size ' num2str(size(parentSubs(1).edges,1)+1) '..']);
        childSubArr = cell(numel(parentSubs),1);
        setDistributions = ones(floor(numel(parentSubs)/numberOfThreads),1) * numberOfThreads;
        if rem(numel(parentSubs), numberOfThreads) > 0
            setDistributions = [setDistributions;rem(numel(parentSubs), numberOfThreads)]; %#ok<AGROW>
        end
        parentSubSets = cell(mat2cell(1:numel(parentSubs), 1, setDistributions));
        for setItr = 1:numel(parentSubSets)
            %% Step 2.1: If it has been too long, we need to finish execution.
            elapsedTime = toc(startTime);
            haltedParent = -1;
            if elapsedTime > maxTime
                haltedParent = parentSubSets{setItr}(1);
                display(['[SUBDUE] Time is up! Breaking from process before ' ...
                    'extending parent ' num2str(haltedParent) ' of size ' num2str(size(parentSubs(1).edges,1)+1) '..']);
                break;
            end
            
            % All good, continue.
            processedSet = parentSubSets{setItr};
            parfor parentItr = processedSet
                %% Step 2.2: Extend head in all possible directions into childSubs.
                childSubs = extendSub(parentSubs(parentItr), allEdges);
                if isempty(childSubs) 
                    continue;
                end
                %% Step 2.3: Remove duplicates from childSubs.
                % Here, we check if childSubs has duplicates in extendedSubs.
    %            childSubs = removeDuplicateSubs([extendedSubs, childSubs]);
                childSubs = childSubs([childSubs.mdlScore] == 0);

                %% Step 2.4: Evaluate childSubs, find their instances.
                childSubs = evaluateSubs(childSubs, evalMetric, allEdges, allSigns, graphSize);

                %% Save childSubs
                childSubArr(parentItr) = {childSubs};
            end
        end
        
        %% Add each children group in childGroupArr into extendedSubs.
        display('[SUBDUE] Merging all children and putting them into bestSubs if they match final criteria.');
        extendedSubs = [];
        for childGroupItr = 1:numel(childSubArr)
            %% Step 2.5: Add childSubs to extendedSubs and bestSubs.
            childSubs = childSubArr{childGroupItr};
            if ~isempty(childSubs)
                extendedSubs = addToQueue(childSubs, extendedSubs, beam);
                if numel(childSubs(1).edges) >= (options.subdue.minSize-1)
                    % Here, a check ensures that the subs to put in bestSubs
                    % are distinct. Simply, an instance of any sub in bestSubs
                    % should not have the exact same node list (children) as
                    % any other sub in the list.
                    [~,sortedIdx]=sort([childSubs.mdlScore]);
                    childSubs=childSubs(sortedIdx);
%                    childSubs = removeEncodedSubs(childSubs, bestSubs, oppositeModes);

                    % Add remaining child subs to best subs.
                    bestSubs = addToQueue(childSubs, bestSubs, nsubs);    
                end 
            end
        end
        
        %% Step 2.6: Swap expandedSubs with parentSubs.
        if haltedParent == -1
            display('[SUBDUE] Swapping children with parents and going on..');
            parentSubs = extendedSubs;
        else
            break;
        end
    end
    
    display('[SUBDUE] Processing has finished. Writing all to output and quitting.');
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
       
       % Assign mdl score.
       vocabLevel(bestSubItr).mdlScore = bestSubs(bestSubItr).mdlScore;
       vocabLevel(bestSubItr).normMdlScore = bestSubs(bestSubItr).normMdlScore;
       
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
       edges = allEdges(centerIdx);
       edgeIdx = {instances.edges}';
       instanceEdges = cellfun(@(x,y) x(y,:), edges, edgeIdx, 'UniformOutput', false);
       instanceChildren = cellfun(@(x,y) [x, y(:,2)'], centerIdxCellArr, instanceEdges, 'UniformOutput', false);
       instanceSigns = num2cell(allSigns(centerIdx));
       
       % Assign fields to graphs.
       [graphLevel(instanceOffset:instanceEndOffset).labelId] = deal(labelIds{:});
       [graphLevel(instanceOffset:instanceEndOffset).children] = deal(instanceChildren{:});
       [graphLevel(instanceOffset:instanceEndOffset).sign] = deal(instanceSigns{:});
       instanceOffset = instanceOffset + numberOfInstances;
       
       % Clear variables.
       clear centerIdx centerIdxCellArr edges instanceEdges instanceChildren instanceSigns
    end
    
   nextVocabLevel = vocabLevel;
   nextGraphLevel = graphLevel;
end


%> Name: getSingleNodeSubs
%>
%> Description: getSingleNodeSubs is used to obtain single-node subs from
%> labels of all instances in allLabels. The result is a number of
%> subs representing compositions each having their instances.
%> 
%> @param allLabels Labels for every graph node.
%> @param allSigns Signs for every graph node.
%>
%> @retval singleNodeSubs Substructure list of single-node subs.
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 24.02.2014
%> Ver 1.1 on 01.09.2014 Removal of global parameters.
function singleNodeSubs = getSingleNodeSubs(allLabels, allSigns)
    numberOfSubs = max(allLabels);
    singleNodeSubs(numberOfSubs) = Substructure();
    validSubs = ones(numberOfSubs,1)>0;
    
    %% For each center node label type, we create a substructure.
    for subItr = 1:numberOfSubs
        subCenterIdx = allLabels == subItr;
        instances = find(subCenterIdx);
        numberOfInstances = numel(instances);
        
        %Assign center id.
        singleNodeSubs(subItr).centerId = subItr;
        
        % Give maximum score so that it is at the top of things.
        singleNodeSubs(subItr).mdlScore = numberOfInstances;
        if numberOfInstances>0
            singleNodeInstances(numberOfInstances) = Instance(); %#ok<AGROW>

            % Fill in instance information. 
            instanceIdx = allLabels == subItr;
            subNodeAssgnArr = num2cell([find(instanceIdx), allSigns(instanceIdx,1)]);
            [singleNodeInstances.centerIdx, singleNodeInstances.sign] = deal(subNodeAssgnArr{:});

            singleNodeSubs(subItr).instances = singleNodeInstances;
            clear singleNodeInstances;
        else
            validSubs(subItr) = 0;
        end
    end
    
    % Eliminate those with no instances.
    singleNodeSubs = singleNodeSubs(validSubs);
end

%> Name: extendSub
%>
%> Description: extendSub(..) extends 'sub' in all possible ways
%> by extending its instances, and returns the new sub list 'extendedSubs'
%> along with instances of each sub in returned list.
%> 
%> @param sub Sub to be extended.
%> @param allEdges List of all edges in the graph.
%>
%> @retval extendedSubs Extended sub list.
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 24.02.2014
%> Ver 1.1 on 01.09.2014 Removal of global parameters.
function [extendedSubs] = extendSub(sub, allEdges)
    centerIdx = cat(1,sub.instances.centerIdx);
    subAllEdges = allEdges(centerIdx);
    % Get unused edges from sub's instances.
    allUsedEdgeIdx = {sub.instances.edges}';
    allUnusedEdges = cellfun(@(x,y) x(setdiff(1:size(x,1),y),:), subAllEdges, allUsedEdgeIdx, 'UniformOutput', false);
    
    % Record which edge belongs to which instance. 
    allEdgeInstanceIds = zeros(size(allUnusedEdges,1),1);
    itrOffset = 1;
    unusedEdgeCount = cellfun(@(x) size(x,1), allUnusedEdges);
    for itr = 1:numel(allUnusedEdges)
        beginOffset = itrOffset;
        allEdgeInstanceIds(beginOffset:(beginOffset+(unusedEdgeCount(itr)-1))) = itr;
        itrOffset = itrOffset + unusedEdgeCount(itr);
    end         
    allUnusedEdges = cat(1, allUnusedEdges{:});
    
    % If no edges remain, exit.
    if isempty(allUnusedEdges)
        extendedSubs = [];
        return; 
    end
    
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
        edgeInstanceIds = allEdgeInstanceIds(edgesToExtendIdx)';
        newInstances = sub.instances(edgeInstanceIds);
        
        % Each added edge actually means a new instance. 
        for instanceItr = 1:numel(edgeInstanceIds)
            instanceEdges = subAllEdges{edgeInstanceIds(instanceItr)};
            if isempty(instanceEdges)
               continue; 
            end
            addedEdgeIdx = find(instanceEdges(:, 2) == edgesToExtend(instanceItr,2));
            newInstances(instanceItr).edges = [newInstances(instanceItr).edges; addedEdgeIdx];
        end
        
        % Assign instances and we are good to go.
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
%> @param evalMetric Evaluation metric.
%> @param allEdges List of all edges in the graph.
%> @param allEdges List of all signs of the nodes in the graph.
%> @param allEdges Size of the graph.
%>
%> @retval subs Evaluated sub list.
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 24.02.2014
%> Ver 1.1 on 01.09.2014 Removal of global parameters.
function [subs] = evaluateSubs(subs, evalMetric, allEdges, allSigns, graphSize)
    numberOfSubs = numel(subs);
    if strcmp(evalMetric, 'mdl')
        
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
            instanceEdges = allEdges(centerIdx);
            instanceSigns = allSigns(cat(1, subs(subItr).instances.centerIdx));
            instanceUsedEdgeIdx = {subs(subItr).instances.edges}';
            
            % Calculate direct neighbors and neighbors of neighbors for mdl
            % calculation.
            instanceAllNeighbors = cellfun(@(x) x(:,2), instanceEdges, 'UniformOutput', false);
            instanceNeighbors = cellfun(@(x,y) x(y,2), instanceEdges, instanceUsedEdgeIdx, 'UniformOutput', false);
            instanceNeighborEdges = cellfun(@(x) cat(1, allEdges{x}), instanceNeighbors, 'UniformOutput', false);
            nonemptyNeighborIdx = ~cellfun(@(x) isempty(x), instanceNeighborEdges);
            instanceSecondaryNeighbors = cell(size(nonemptyNeighborIdx,1),1);
            instanceSecondaryNeighbors(nonemptyNeighborIdx) = cellfun(@(x) x(:,2), instanceNeighborEdges(nonemptyNeighborIdx), 'UniformOutput', false);
            
            % Learn number of instances/edges in each instance.
            numberOfInstances = size(instanceSigns,1);
            numberOfEdgesInSub = size(subs(subItr).edges,1);
            positiveInstanceIdx = instanceSigns == 1;
            negativeInstanceIdx = instanceSigns == 0;
            numberOfPositiveInstances = numel(find(positiveInstanceIdx));
            
            % 1) Calculate edge-based DL reduction for each instance. 
            % (numel(fastintersect(x,y)) : Number of edges to be removed to
            % represent edges from neighbors to secondary neighbors, since 
            % we connect them directly to the center node.
            % (numel(y) - numel(unique(y))) : Number of duplicates within
            % y. Each secondary neighbor (neighbor of neighbor) needs to be
            % represented only once. 
            edgeDLReductions = cellfun(@(x,y) numel(fastintersect(x,y)) + (numel(y) - numel(unique(y))), instanceAllNeighbors, instanceSecondaryNeighbors);
            
            % 2) Calculate edge-based DL reduction amount due to the removal
            % of edges with each instance. 
            % amount = (nPos - (nNeg - nPos)) * edgePerSub = (2 * nPos - nNeg) * edgePerSub.
            selfEdgeDLReductionAmount = (2*numberOfPositiveInstances-numberOfInstances) * numberOfEdgesInSub;
            
            % 3) Calculate # of nodes to be removed from the graph.
            selfNodeDLReductionAmount = numel(unique([instanceNeighbors{positiveInstanceIdx}])) - ...
                numel(unique([instanceNeighbors{negativeInstanceIdx}]));
            
            %% Estimate DL reduction. (CHECK AGAIN FOR CORRECTNESS!)
            dlDiff = 2*(sum(edgeDLReductions(positiveInstanceIdx)) - sum(edgeDLReductions(negativeInstanceIdx)) + ...
                selfEdgeDLReductionAmount) + selfNodeDLReductionAmount - (numberOfEdgesInSub * 3);
            subs(subItr).mdlScore = dlDiff;
            subs(subItr).normMdlScore = 1 - (dlDiff / graphSize);
        end
    else
        for subItr = 1:numberOfSubs
            subs(subItr).mdlScore = size(subs(subItr).edges,1) * (numel(subs(subItr).instances)-1);
        end
    end
end

% %> Name: removeDuplicateSubs
% %>
% %> Description: Given a substructure list, this function removes duplicates
% %> from it by calculating a unique description for each sub, and then
% %> checking for redundant descriptions (a.k.a. duplicate subs).
% %> 
% %> @param childSubs The substructure list from which duplicates will be
% %> removed.
% %>
% %> @retval childSubs Unique substructures.
% %> 
% %> Author: Rusen
% %>
% %> Updates
% %> Ver 1.0 on 24.02.2014
% function childSubs = removeDuplicateSubs(childSubs)
%     if isempty(childSubs)
%        return; 
%     end
%     edges = {childSubs.edges};
%     edges = cellfun(@(x) sortrows(x), edges, 'UniformOutput', false);
%     edges = cellfun(@(x) (x(:)), edges, 'UniformOutput', false);
%     edges = cell2mat(edges)';
%     subIdentifiers = [cat(1, childSubs.centerId), edges];
%     [~, uniqueSubIdentifiers, ~] = unique(subIdentifiers, 'rows', 'stable');
%     childSubs = childSubs(uniqueSubIdentifiers);
% end

% %> Name: removeEncodedSubs
% %>
% %> Description: Given subs to check (subs2Check) and already considered
% %> subs (encodedSubs), this function removes the subs from subs2Check and
% %> returns remaining subs in validSubs.
% %> 
% %> @param subs2Check Subs each of which is considered for elimination if
% %> they already exist either in encodedSubs, or subs2Check as duplicates.
% %> @param encodedSubs Subs already encoded.
% %> @param oppositeEdgeLabelList nx1 array indicating an edge type's reverse
% %> geometric information if it is encoded by another edge. 
% %> oppositeEdgeLabels(x)=y => oppositeEdgeLabels(y)=x.
% %>
% %> @retval validSubs Unique subs from subs2Check.
% %> 
% %> Author: Rusen
% %>
% %> Updates
% %> Ver 1.0 on 15.02.2014
% function validSubs = removeEncodedSubs(subs2Check, encodedSubs, oppositeEdgeLabelList)
%     if isempty(subs2Check) || isempty(oppositeEdgeLabelList) 
%         validSubs = subs2Check;
%         return;
%     end
%     %% First, eliminate duplicate subs in subs2Check.
%     % Check number of nodes in subs2Check. It is assumed that all subs in
%     % subs2Check have the same size. If size ~= 2, there is no ambiguity in
%     % the geometry, and each sub encodes a different structure. No need to
%     % check.
%     numberOfNodes = size(subs2Check(1).edges,1) + 1;
%     
%     if numberOfNodes == 2
%         % First, we check for duplicate subs in subs2Check.
%         centerIds = cat(1, subs2Check.centerId);
%         edges = {subs2Check.edges};
%         descriptions = [centerIds, cat(1, edges{:})];
%         reverseDescriptions = [descriptions(:,3), oppositeEdgeLabelList(descriptions(:,2)), descriptions(:,1)];
% 
%         % Get description matches.
%         [~, reverseIdx] = ismember(reverseDescriptions, descriptions, 'rows');
% 
%         % Eliminate half of the matches, since one of the two matching subs
%         % will be left in the sub list.
%         reverseCheckArr = (1:numel(reverseIdx))';
%         remainingSubIdx = reverseIdx <= reverseCheckArr;
%         subs2Check = subs2Check(remainingSubIdx);
%     else
%         remainingSubIdx = ones(numel(subs2Check),1) > 0;
%     end
%     validSubs = subs2Check;
%     
%     if isempty(encodedSubs)
%         return;
%     end
%     %% Second task involves checking subs2Check against encodedSubs.
%     encodedSubsNodeCountArr = {encodedSubs.edges};
%     encodedSubsNodeCountArr = cellfun(@(x) size(x,1), encodedSubsNodeCountArr)' + 1;
%     encodedSubs = encodedSubs(encodedSubsNodeCountArr == numberOfNodes);
%     
%     if numberOfNodes == 2
%         % Eliminated those not having exact same number of nodes. If empty,
%         % return.
%         if isempty(encodedSubs)
%             return
%         end
%         
%         % Get reverse descriptions of each sub in subs2Check.
%         reverseDescriptions = reverseDescriptions(remainingSubIdx, :);
% 
%         % Get descriptions of each sub in encodedSubs.
%         encodedCenterIds = cat(1, encodedSubs.centerId);
%         encodedEdges = {encodedSubs.edges};
%         encodedDescriptions = [encodedCenterIds, cat(1, encodedEdges{:})];
% 
%         % Eliminate those in subs2Check that match already encoded structure in
%         % encodedSubs.
%         remainingSubIdx = ~ismember(reverseDescriptions, encodedDescriptions, 'rows');
%         subs2Check = subs2Check(remainingSubIdx);
%     end
%     validSubs = subs2Check;
%     
% end

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


