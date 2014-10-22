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
%> Ver 1.2 on 02.09.2014 Adding display commentary.
function [nextVocabLevel, nextGraphLevel, prevGraphData] = runSubdue(vocabLevel, graphLevel, options)
    %% First thing we do is to convert vocabLevel and graphLevel into different data structures.
    % This process is done to assure fast, vectorized operations.
    % Initialize the priority queue.
    
    bestSubs = [];
    parentSubs = [];
    nextVocabLevel = [];
    nextGraphLevel = [];
    prevGraphData = [];
    
    %% Get the parameters.
    evalMetric = options.subdue.evalMetric;
    mdlNodeWeight = options.subdue.mdlNodeWeight;
    mdlEdgeWeight = options.subdue.mdlEdgeWeight;
    isMDLExact = options.subdue.isMDLExact;
    overlap = options.subdue.overlap;
    beam = options.subdue.beam;
    nsubs = options.subdue.nsubs;
    maxTime = options.subdue.maxTime;
    maxSize = options.subdue.maxSize;
    numberOfThreads = options.numberOfThreads;
    
    %% Initialize data structures.
    display('[SUBDUE] Initializing data structures for internal use..');
    % Helper data structures.
    assignedEdges = {graphLevel.adjInfo};    
    if isempty(assignedEdges) 
        return;
    else
        allEdgeNodePairs = cat(1,assignedEdges{:});
        if isempty(allEdgeNodePairs)
            return;
        end
        numberOfAllEdges = size(allEdgeNodePairs,1);
        clear allEdgeNodePairs;
    end
    nonemptyEdgeIdx = cellfun(@(x) ~isempty(x), assignedEdges);
    allLabels = cat(1, graphLevel.labelId);
    assignedEdges(nonemptyEdgeIdx) = cellfun(@(x) [x(:,1:3), allLabels(x(:,2))], ...
        assignedEdges(nonemptyEdgeIdx), 'UniformOutput', false);
    allEdges(numel(graphLevel)) = AdjInfo();
    [allEdges.adjInfo] = deal(assignedEdges{:});
    allEdgeNodePairs = cat(1,allEdges.adjInfo);
    allEdgeNodePairs = allEdgeNodePairs(:,1:2);
    clear assignedEdges;
    
    % If no edges are present, time to return.
    allSigns = cat(1, graphLevel.sign);
    
    % Graph size formulation is very simple: edgeWeight * #edges + edgeWeight * #nodes. 
    graphSize = numberOfAllEdges * mdlEdgeWeight + ...
        numel(graphLevel) * mdlNodeWeight;
    
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

                %% Step 2.3: Evaluate childSubs, find their instances.
                childSubs = evaluateSubs(childSubs, evalMetric, allEdges, allEdgeNodePairs, ...
                    allSigns, graphSize, overlap, mdlNodeWeight, mdlEdgeWeight, isMDLExact);

                %% Save childSubs
                childSubArr(parentItr) = {childSubs};
            end
        end
        
        %% Add each children group in childGroupArr into extendedSubs.
        display('[SUBDUE] Merging all children and putting them into bestSubs if they match final criteria.');
        extendedSubs = [];
        %% Step 2.4: Add childSubs to extendedSubs and bestSubs.
        childSubArr = cat(2,childSubArr{:});
        
        % Remove duplicate nodes from the children subs.
        if size(childSubArr(1).edges,1) > 1
            childSubEdges = cell(1, numel(childSubArr));
            for subItr = 1:numel(childSubArr)
                childSubEdges(subItr) = {sortrows(childSubArr(subItr).edges)};
            end
            subCenters = {childSubArr.centerId};
            vocabDescriptors = cellfun(@(x,y) num2str([x; y(:)]'), subCenters, childSubEdges, 'UniformOutput', false);
            [~, validChildrenIdx, ~] = unique(vocabDescriptors, 'stable');
            childSubArr = childSubArr(validChildrenIdx);
        end
        
        % Add children to the both queues.
        if ~isempty(childSubArr)
            extendedSubs = addToQueue(childSubArr, extendedSubs, beam);
            bestSubs = addToQueue(childSubArr, bestSubs, nsubs); 
        end
        clear childSubArr;
        
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
    
    %% Previous graph data is returned.
    prevGraphData.allEdges = allEdges;
    prevGraphData.allEdgeNodePairs = allEdgeNodePairs;
    prevGraphData.allSigns = allSigns;
    prevGraphData.bestSubs = bestSubs;
    
    %% Allocate space for new graphLevel and vocabLevel.
    numberOfInstances = 0;
    for bestSubItr = 1:numberOfBestSubs
        numberOfInstances = numberOfInstances + numel(bestSubs(bestSubItr).instances);
    end
    clear vocabLevel graphLevel
    if numberOfInstances>0
        vocabLevel(numberOfBestSubs) = options.vocabNode;
        graphLevel(numberOfInstances) = options.graphNode;

        %% Assign instances of each discovered substructure to graph node. 
        % We store substructures in vocabLevel array.
        instanceOffset = 1;
        for bestSubItr = 1:numberOfBestSubs
           % Assign label of sub.
           vocabLevel(bestSubItr).label = int32(bestSubItr);

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
           labelIds = num2cell(repmat(int32(bestSubItr), numberOfInstances,1));

           %% Create required fields such as centerIdx, edges, children and sign array to assign to instances.
           centerIdx = [instances.centerIdx];
           centerIdxCellArr = num2cell(centerIdx');
           numberOfInstances = numel(centerIdx);
           instanceEndOffset = instanceOffset + numberOfInstances - 1;
           edges = {allEdges(centerIdx).adjInfo}';
           edgeIdx = {instances.edges}';
           instanceEdges = cellfun(@(x,y) x(y,:), edges, edgeIdx, 'UniformOutput', false);
           instanceChildren = cellfun(@(x,y) [x, y(:,2)'], centerIdxCellArr, instanceEdges, 'UniformOutput', false);
           instanceSigns = num2cell(allSigns(centerIdx));

           % Assign fields to graphs.
           [graphLevel(instanceOffset:instanceEndOffset).labelId] = deal(labelIds{:});
           [graphLevel(instanceOffset:instanceEndOffset).children] = deal(instanceChildren{:});
           [graphLevel(instanceOffset:instanceEndOffset).sign] = deal(instanceSigns{:});
           instanceOffset = instanceOffset + numberOfInstances;
        end
    else
        vocabLevel = [];
        graphLevel = [];
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
            subNodeAssgnArr = num2cell([int32(find(instanceIdx)), allSigns(instanceIdx,1)]);
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
%> @param allEdges List of all edges in the graph (Nx4).
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
    subAllEdges = {allEdges(centerIdx).adjInfo}';
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
%> Key points:
%>  1) DL estimation not rigorous. It is considered that 
%>        each node needs a label (int) and a pointer to its edges (int)
%>        each edge needs a node label (int) for its destination node, an
%>        edge label(int) and a binary isDirected label (bit)
%>     in the final graph description. Node and edge weights in DL
%>     calculation is stored in options.
%> 
%> @param subs Sub list which will be evaluated.
%> @param evalMetric Evaluation metric.
%> @param allEdges List of all edges in the graph.
%> @param allEdgeNodePairs List of all edge node pairs in the graph.
%> @param allSigns List of all signs of the nodes in the graph.
%> @param graphSize Size of the graph.
%> @param overlap If true, overlapping instances are considered in
%> evaluation of a sub.
%> @param mdlNodeWeight Node weight in DL calculations (MDL).
%> @param mdlEdgeWeight Edge weight in DL calculations (MDL).
%> @param isMDLExact If true, exact MDL calculation. Approximate otherwise.
%>
%> @retval subs Evaluated sub list.
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 24.02.2014
%> Ver 1.1 on 01.09.2014 Removal of global parameters.
function [subs] = evaluateSubs(subs, evalMetric, allEdges, allEdgeNodePairs, allSigns, graphSize, overlap, mdlNodeWeight, mdlEdgeWeight, isMDLExact)
    numberOfSubs = numel(subs);
    for subItr = 1:numberOfSubs
        % We compress the object graph using the children, and the
        % edges they are involved. 
        subScore = getSubScore(subs(subItr), allEdges, allEdgeNodePairs, evalMetric, ...
           allSigns, mdlNodeWeight, mdlEdgeWeight, ....
            overlap, isMDLExact);

        %% Assign the score of the sub, as well as its normalized mdl score if applicable.
        subs(subItr).mdlScore = subScore;
        if strcmp(evalMetric, 'mdl')
            subs(subItr).normMdlScore = 1 - (subScore / graphSize);
        end
    end
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