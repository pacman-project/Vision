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
function [nextVocabLevel, nextGraphLevel, optimalThreshold] = runSubdue(vocabLevel, graphLevel, nodeDistanceMatrix, edgeDistanceMatrix, categoryArrIdx, validationIdx, options)
    %% First thing we do is to convert vocabLevel and graphLevel into different data structures.
    % This process is done to assure fast, vectorized operations.
    % Initialize the priority queue.
    
    bestSubs = [];
    parentSubs = [];
    nextVocabLevel = [];
    nextGraphLevel = [];
    
    %% Get the parameters.
    evalMetric = options.subdue.evalMetric;
    mdlNodeWeight = options.subdue.mdlNodeWeight;
    mdlEdgeWeight = options.subdue.mdlEdgeWeight;
    isMDLExact = options.subdue.isMDLExact;
    overlap = options.subdue.overlap;
    beam = options.subdue.beam;
    nsubs = options.subdue.nsubs;
    maxTime = options.subdue.maxTime;
    minSize = options.subdue.minSize;
    maxSize = options.subdue.maxSize;
    isSupervised = options.subdue.supervised;
    orgThreshold = single(options.subdue.threshold);
    maxThreshold = options.subdue.maxThreshold;
    minThreshold = options.subdue.minThreshold;
    maxDepth = options.subdue.thresholdSearchMaxDepth;
    singlePrecision = options.singlePrecision;
    reconstructionFlag = options.reconstruction.flag;
    stoppingCoverage = options.reconstruction.stoppingCoverage;
    optimalThreshold = orgThreshold;
    noveltyThr = options.noveltyThr;
    if options.validationFlag
        validationFolds = options.validationFolds;
    else
        validationFolds = 1; 
    end
    if reconstructionFlag
        singleNodeThreshold = maxThreshold + singlePrecision;
    else
        singleNodeThreshold = orgThreshold +singlePrecision; % Hard threshold for cost of matching two subs.
    end
    
    % At this point we get more subs than we need, since we're trying to
    % optimize based on the number of subs.
    numberOfReconstructiveSubs = options.reconstruction.numberOfReconstructiveSubs;
    
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
    prevActivations = cat(1, graphLevel.activation);
    assignedEdges(nonemptyEdgeIdx) = cellfun(@(x) [x(:,1:3), allLabels(x(:,2))], ...
        assignedEdges(nonemptyEdgeIdx), 'UniformOutput', false);
    allEdges(numel(graphLevel)) = AdjInfo();
    [allEdges.adjInfo] = deal(assignedEdges{:});
    allEdgeNodePairs = cat(1,allEdges.adjInfo);
    allEdgeNodePairs = allEdgeNodePairs(:,1:2);
    clear assignedEdges;
    
    % If no edges are present, time to return.
    allSigns = uint8(cat(1, graphLevel.sign));
    categoryArrIdx = uint8(categoryArrIdx(cat(1, graphLevel.imageId)))';
    validationIdx = validationIdx(cat(1, graphLevel.imageId));
    
    % Graph size formulation is very simple: edgeWeight * #edges + edgeWeight * #nodes. 
    graphSize = numberOfAllEdges * mdlEdgeWeight + ...
        numel(graphLevel) * mdlNodeWeight;
    
%    % Find the total cost of matching considering the max size, and set the
%    % threshold.
   if reconstructionFlag
        adaptiveThreshold = maxThreshold * (single(maxSize) * 2 - 1) + singlePrecision;
   else
        adaptiveThreshold = orgThreshold * (single(maxSize) * 2 - 1) + singlePrecision;
   end
   
    %% Step 1:Find single-vertex subs, and put them into beamSubs.
    display('[SUBDUE] Creating single node subs..');
    singleNodeSubs = getSingleNodeSubs(allLabels, allSigns, ...
        nodeDistanceMatrix, categoryArrIdx, validationIdx, adaptiveThreshold);
    parentSubs = addToQueue(singleNodeSubs, parentSubs, numel(singleNodeSubs));
    singleNodeSubsFinal = [];
    
    %% Step 2: Main loop
    startTime = tic;
    % Let's find the adaptive threshold with the size of the minimum 
    currentSize = 2;
    if reconstructionFlag
        overallThreshold = maxThreshold * (single(maxSize) * 2 - 1) + singlePrecision;
    else
        overallThreshold = orgThreshold * (single(maxSize) * 2 - 1) + singlePrecision;
    end
    minMdlScoreFinal = -inf; % This two will increase as more and more subs are discovered.
    bestMdlScoresFinal = [];
    
    while ~isempty(parentSubs)
        minMdlScoreExtend = -inf;   
        % This two will increase as more and more subs are discovered.
        % If we're building a reconstructive hierarchy, we do not want to
        % eliminate subs, as we may need them later. That's why we are
        % setting adaptiveMinThreshold to the minimum value.
        % adaptiveThreshold is used to eliminate the subs to be put in the
        % final node list.
        if reconstructionFlag
            adaptiveThreshold = maxThreshold * (single(currentSize) * 2 - 1) + singlePrecision;
            adaptiveMinThreshold = minThreshold * (single(currentSize) * 2 - 1) + singlePrecision;
        else
            adaptiveThreshold = orgThreshold * (single(currentSize) * 2 - 1) + singlePrecision;
            adaptiveMinThreshold = adaptiveThreshold;
        end
        
        if ~isempty(bestSubs)
            bestMdlScoresFinal = [bestSubs.mdlScore];
        end
        
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
        
        %% Evaluate single node subs, if required.
        if minSize == 1
            if currentSize == 2
                 % Find single node subs.
                 singleNodeSubsFinal = getSingleNodeSubs(allLabels, allSigns, ...
                    nodeDistanceMatrix, categoryArrIdx, validationIdx, singleNodeThreshold);
               
               % Evaluate them.
                singleNodeSubsFinal = evaluateSubs(singleNodeSubsFinal, evalMetric, allEdges, allEdgeNodePairs, ...
                    allSigns, graphSize, overlap, mdlNodeWeight, mdlEdgeWeight, isMDLExact, isSupervised);

                %% Remove those with no instances. 
                centerIdxArr = {singleNodeSubsFinal.instanceCenterIdx};
                validSingleSubs = cellfun(@(x) ~isempty(x), centerIdxArr);
                singleNodeSubsFinal = singleNodeSubsFinal(validSingleSubs);

                %% Sort singleNodeSubsFinal by their scores.
                % We do not eliminate nodes with negative mdl scores,
                % as single node subs do not provide any compression.
                mdlScores = [singleNodeSubsFinal.mdlScore];
                [~, sortIdx] = sort(mdlScores, 'descend');
                singleNodeSubsFinal = singleNodeSubsFinal(sortIdx);
            end
        else
            singleNodeSubsFinal = [];
        end
        
        %% All good, execution. First, we put parents into fixed size subsets for easy parallelization.
        display(['[SUBDUE/Parallel] Expanding subs of size ' num2str(size(parentSubs(1).edges,1)+1) '..']);
        childSubArrFinal = cell(numel(parentSubs),1);
        childSubArrExtend = cell(numel(parentSubs),1);
        mdlScoreArrFinal = cell(numel(parentSubs),1);
        mdlScoreArrExtend = cell(numel(parentSubs),1);
        numberOfParentSubs = numel(parentSubs);
        if numberOfParentSubs > 50
            setDistributions = 1:50:numberOfParentSubs;
            if setDistributions(end) ~= numberOfParentSubs
                setDistributions = [setDistributions, numberOfParentSubs]; %#ok<AGROW>
            end
            setDistributions = setDistributions(2:(end)) - setDistributions(1:(end-1));
            setDistributions(end) = setDistributions(end) + 1;
        else
            setDistributions = numberOfParentSubs;
        end
        parentSubSets = cell(mat2cell(1:numel(parentSubs), 1, setDistributions));
        for setItr = 1:numel(parentSubSets)
            %% Step 2.1: If it has been too long, we need to finish execution.
            elapsedTime = toc(startTime);
            haltedParent = -1;
            if elapsedTime > maxTime && currentSize > 2
                haltedParent = parentSubSets{setItr}(1);
                display(['[SUBDUE] Time is up! Breaking from process before ' ...
                    'extending parent ' num2str(haltedParent) ' of size ' num2str(size(parentSubs(1).edges,1)+1) '..']);
                break;
            end
            
            %% All good, continue with the main algorithm.
            processedSet = parentSubSets{setItr};
            parfor parentItr = processedSet
                %% Step 2.2: Extend head in all possible directions into childSubs.
                display(['[SUBDUE/Parallel] Expanding sub ' num2str(parentItr) ' of size ' num2str(currentSize-1) '..']);
                childSubs = extendSub(parentSubs(parentItr), allEdges, nodeDistanceMatrix, edgeDistanceMatrix, ...
                    singlePrecision, overallThreshold);
                if isempty(childSubs) 
                    continue;
                end
                
                %% Step 2.3: Eliminate subs that match to more frequent subs.
                [childSubs, ~] = getNonOverlappingSubs(childSubs, nodeDistanceMatrix, edgeDistanceMatrix, ...
                     adaptiveMinThreshold, singlePrecision);

                 % Get the list of final subs and indices of subs chosen
                 % for extension. 
                [childSubsFinal, childSubsExtend] = getFinalSubs(childSubs, adaptiveThreshold);
                 
                %% Step 2.4: Evaluate childSubs, find their instances.
                childSubsFinal = evaluateSubs(childSubsFinal, evalMetric, allEdges, allEdgeNodePairs, ...
                    allSigns, graphSize, overlap, mdlNodeWeight, mdlEdgeWeight, isMDLExact, isSupervised);
                
                % Assign mdl scores of subs chosen for extension as well. 
                [childSubsFinal, childSubsExtend] = copyMdlScores(childSubsFinal, childSubsExtend);
                
                %% Eliminate childSubs which will render useless to preserve memory.
                mdlScores = [childSubsFinal.mdlScore];
                validMdlScoreIdxFinal = mdlScores > minMdlScoreFinal;
                childSubsFinal = childSubsFinal(validMdlScoreIdxFinal);
                validMdlScoreIdxExtend = mdlScores > minMdlScoreExtend;
                childSubsExtend = childSubsExtend(validMdlScoreIdxExtend);
                    
                % If no subs are remaining, continue.
                if isempty(childSubsFinal) 
                    continue;
                end
                
                %% Sort childSubs by the mdl scores.
                mdlScores = [childSubsFinal.mdlScore];
                [sortedMdlScoresFinal, sortIdx] = sort(mdlScores, 'descend');
                childSubsFinal = childSubsFinal(sortIdx);
                
                mdlScores = [childSubsExtend.mdlScore];
                [sortedMdlScoresExtend, sortIdx] = sort(mdlScores, 'descend');
                childSubsExtend = childSubsExtend(sortIdx);
                
                %% Save childSubs and extended subs.
                childSubArrFinal(parentItr) = {childSubsFinal};
                childSubArrExtend(parentItr) = {childSubsExtend};
                mdlScoreArrFinal(parentItr) = {sortedMdlScoresFinal};
                mdlScoreArrExtend(parentItr) = {sortedMdlScoresExtend};
            end
            
            %% Based on the best mdl scores, update minMdlScoreFinal and minMdlScoreExtend.
            newMdlScores = [bestMdlScoresFinal, [mdlScoreArrFinal{:}]];
            if numel(newMdlScores) > nsubs
                newMdlScores = sort(newMdlScores, 'descend');
                minMdlScoreFinal = newMdlScores(nsubs);
            end
            newMdlScores = [mdlScoreArrExtend{:}];
            if numel(newMdlScores) > beam
                newMdlScores = sort(newMdlScores, 'descend');
                minMdlScoreExtend = newMdlScores(beam);
            end
            
            % In addition, remove the subs that have low mdl scores.
            % First, we handle final subs.
            if minMdlScoreFinal ~= -inf
                validBestSubIdx = bestMdlScoresFinal >= minMdlScoreFinal;
                bestSubs = bestSubs(validBestSubIdx);
                nonemptyArrIdx = cellfun(@(x) ~isempty(x), childSubArrFinal);
                childSubArrFinal(nonemptyArrIdx) = cellfun(@(x, y) x(y >= minMdlScoreFinal), childSubArrFinal(nonemptyArrIdx), mdlScoreArrFinal(nonemptyArrIdx), 'UniformOutput', false);
                mdlScoreArrFinal(nonemptyArrIdx) = cellfun(@(x) x(x >= minMdlScoreFinal), mdlScoreArrFinal(nonemptyArrIdx), 'UniformOutput', false);
            end
            
            % Then, we handle extended subs.
            if minMdlScoreExtend ~= -inf
                nonemptyArrIdx = cellfun(@(x) ~isempty(x), childSubArrExtend);
                childSubArrExtend(nonemptyArrIdx) = cellfun(@(x, y) x(y >= minMdlScoreExtend), childSubArrExtend(nonemptyArrIdx), mdlScoreArrExtend(nonemptyArrIdx), 'UniformOutput', false);
                mdlScoreArrExtend(nonemptyArrIdx) = cellfun(@(x) x(x >= minMdlScoreExtend), mdlScoreArrExtend(nonemptyArrIdx), 'UniformOutput', false);
            end
        end
        
        %% Add each children group in childGroupArr into extendedSubs.
        display('[SUBDUE] Merging all children and putting them into bestSubs if they match final criteria.');
        extendedSubs = [];
        %% Step 2.4: Add childSubs to extendedSubs and bestSubs.
        childSubArrFinal = cat(2,childSubArrFinal{:});
        childSubArrExtend = cat(2, childSubArrExtend{:});
        % Add children to both queues.
        if ~isempty(childSubArrFinal)
            %% A substructure has many different parse trees at this point.
            % We remove duplicate subs from childSubArr.
            [childSubArrFinal] = removeDuplicateSubs(childSubArrFinal);
            
            % Check for minSize to put to final sub array.
            if currentSize >= minSize
                bestSubs = addToQueue(childSubArrFinal, bestSubs, nsubs);
                if numel(bestSubs) == nsubs
                    minMdlScoreFinal = min([bestSubs.mdlScore]);
                end
            end
        end
        if ~isempty(childSubArrExtend)
            %% A substructure has many different parse trees at this point.
            % We remove duplicate subs from childSubArr.
            [childSubArrExtend] = removeDuplicateSubs(childSubArrExtend);
            % Check for maxSize to put to extendedSubs (parentSubs for next level).
            if currentSize < maxSize
                extendedSubs = addToQueue(childSubArrExtend, extendedSubs, beam);
            end
        end
            
        clear childSubArrExtend childSubArrFinal;
        
        %% Step 2.6: Swap expandedSubs with parentSubs.
        if haltedParent == -1
            display('[SUBDUE] Swapping children with parents and going on..');
            parentSubs = extendedSubs;
            currentSize = currentSize + 1;
            clear extendedSubs;
        else
            break;
        end
    end
    
    %% At this point, we add single node subs to the list of known subs.
    if ~isempty(singleNodeSubsFinal)
        bestSubs = addToQueue(singleNodeSubsFinal, bestSubs, nsubs + numel(singleNodeSubsFinal)); 
    end
    %% Step 3: Create nextVocabLevel and nextGraphLevel from bestSubs.
    numberOfBestSubs = numel(bestSubs);
    
    % If no subs detected, return.
    if numberOfBestSubs < 1
       return;
    end
    
    % Remove duplicate instances from each sub.
    display('[SUBDUE/Parallel] Removing duplicate instances from each sub, and updating their match scores..');
    bestSubs = removeDuplicateInstances(bestSubs);
    
    %% Allocate space for new graphLevel and vocabLevel.
    numberOfInstances = 0;
    for bestSubItr = 1:numberOfBestSubs
        numberOfInstances = numberOfInstances + size(bestSubs(bestSubItr).instanceCenterIdx,1);
    end
    clear vocabLevel graphLevel
    if numberOfInstances>0
       %% If required, we'll pick best parts based on the reconstruction of the data.
       if reconstructionFlag
           [bestSubs, optimalThreshold] = selectPartsUnsupervised(bestSubs, ...
               nodeDistanceMatrix, edgeDistanceMatrix, singlePrecision, ...
               stoppingCoverage, numberOfReconstructiveSubs, minThreshold, maxThreshold, ...
               maxDepth, validationFolds);
           
           bestSubs = evaluateSubs(bestSubs, 'mdl', allEdges, allEdgeNodePairs, ...
               allSigns, graphSize, overlap, mdlNodeWeight, mdlEdgeWeight, 1, isSupervised); 
%            %% For now, we do mdl-based sorting. The mdl scores are normalized by the size of the bestSub.
%            for bestSubItr = 1:numel(bestSubs)
%                 bestSubs(bestSubItr).mdlScore = bestSubs(bestSubItr).mdlScore / numel(bestSubs(bestSubItr).instanceCenterIdx);
%            end
           
           % Sort bestSubs by their mdl scores.
           mdlScores = [bestSubs.mdlScore];
           [~, mdlSortIdx] = sort(mdlScores, 'descend');
           bestSubs = bestSubs(mdlSortIdx);
       end
       
       %% Learn activations for instances, and save the children for future use.
        display('[SUBDUE/Parallel] Calculating activations for each instance, and saving instance nodes for later use..');
        numberOfBestSubs = numel(bestSubs);
        % Learn the maximum children count in instances.
        maxSubSize = 1;
        for bestSubItr = 1:numberOfBestSubs
            maxSubSize = max(maxSubSize, (size(bestSubs(bestSubItr).edges,1)+1));
        end

        multiplier = round(1/singlePrecision);
        % Find the children for each instance, and use it as a descriptor.
        % We'll then find unique descriptors, which is the initial phase of
        % inhibition.
        instanceChildrenDescriptors = cell(numberOfBestSubs,1);
        instanceActivations = cell(numberOfBestSubs,1);
        remainingInstanceLabels = cell(numberOfBestSubs,1);
        parfor bestSubItr = 1:numberOfBestSubs
           instanceChildren = bestSubs(bestSubItr).instanceChildren;
           tempNum = size(instanceChildren,1);
           % Find children descriptors and save 'em.
           childrenDescriptors = zeros(tempNum, maxSubSize, 'int32');
           childrenDescriptors(:, 1:size(instanceChildren,2)) = instanceChildren;
           instanceChildrenDescriptors{bestSubItr} = childrenDescriptors;
           adaptiveThreshold = single(optimalThreshold * ((size(bestSubs(bestSubItr).edges,1)+1)*2-1)) + singlePrecision;
           instanceMatchScores = fix(multiplier * ((adaptiveThreshold - bestSubs(bestSubItr).instanceMatchCosts) / adaptiveThreshold)) / multiplier;
           
           % Calculate activations for the next level.
           instancePrevActivations = prevActivations(instanceChildren); %#ok<PFBNS>
           if size(instanceChildren,1)>1
               instanceMeanActivations = mean(instancePrevActivations,2);
           else
               instanceMeanActivations = mean(instancePrevActivations);
           end
           activations = instanceMatchScores .* instanceMeanActivations;
           instanceActivations{bestSubItr} = activations;
           
           % Find instance labels.
           instanceLabels = repmat(int32(bestSubItr), tempNum, 1);
           remainingInstanceLabels{bestSubItr} = instanceLabels;
        end
        instanceChildrenDescriptors = cat(1, instanceChildrenDescriptors{:});
        
        %% If inhibition is applied, we eliminate the instances which exist in multiple subs.
        % For each unique instance (node set), we find the best matching
        % sub, and assign the instance to that sub.
        if noveltyThr > 0
           display('[SUBDUE] Performing initial inhibition. Eliminating duplicate instances among different subs..');
           %        Alternative 2
           instanceActivations = cat(1, instanceActivations{:});
           remainingInstanceLabels = cat(1, remainingInstanceLabels{:});
           [~, sortIdx] = sort(instanceActivations, 'descend');
           orderedInstanceChildrenDescriptors = instanceChildrenDescriptors(sortIdx,:);

           % Update data structures by selecting unique children.
           % Find unique rows, which correspond to unique instances.
           [~, IA, ~] = unique(orderedInstanceChildrenDescriptors, 'rows', 'stable');
           IA = sort(sortIdx(IA));
           
           % Get the remaining sub count.
           remainingBestSubs = unique(remainingInstanceLabels(IA))';
           numberOfBestSubs = numel(remainingBestSubs);
           numberOfInstances = numel(IA);
           instanceChildrenDescriptors = instanceChildrenDescriptors(IA,:);
        else
           %        Alternative 1
           remainingBestSubs = 1:numberOfBestSubs;
           IA = 1:size(instanceChildrenDescriptors,1);
           numberOfInstances = numel(IA);
        end
        
        %% Fill in vocabLevel and graphLevel.
       %Allocate space for new graph/vocab level.
       display('[SUBDUE] Processing has finished. Writing all to output and quitting.');
       vocabLevel(numberOfBestSubs) = options.vocabNode;
       graphLevel(numberOfInstances) = options.graphNode;
       
        % First, we start with vocabLevel.
        for bestSubItr = 1:numberOfBestSubs
           % Assign label of sub.
           vocabLevel(bestSubItr).label = int32(bestSubItr);

           % Assign mdl score.
           vocabLevel(bestSubItr).mdlScore = bestSubs(remainingBestSubs(bestSubItr)).mdlScore;
           vocabLevel(bestSubItr).normMdlScore = bestSubs(remainingBestSubs(bestSubItr)).normMdlScore;

           % Assign sub edges (definition)
           numberOfEdges = size(bestSubs(remainingBestSubs(bestSubItr)).edges,1);
           if numberOfEdges > 0
               subEdges = [ones(numberOfEdges,1), ...
                   (2:1:(numberOfEdges+1))', bestSubs(remainingBestSubs(bestSubItr)).edges(:,1), ones(numberOfEdges,1)];
               vocabLevel(bestSubItr).adjInfo = subEdges;
               vocabLevel(bestSubItr).children = [bestSubs(remainingBestSubs(bestSubItr)).centerId; bestSubs(remainingBestSubs(bestSubItr)).edges(:,2)]';
           else
               vocabLevel(bestSubItr).children = bestSubs(remainingBestSubs(bestSubItr)).centerId;
           end
        end
        
        % Now, we fill in the info of vocabLevel.
        actualInstanceOffset = 1;
        instanceOffset = 1;
        remainingSubItr = 1;
        multiplier = round(1/singlePrecision);
        for bestSubItr = 1:numel(bestSubs)
           if ~ismember(bestSubItr,remainingBestSubs)
               instanceOffset = instanceOffset + numel(bestSubs(bestSubItr).instanceCenterIdx);
           else
               %% Assign instances to this sub.
               instanceEndOffset = instanceOffset + numel(bestSubs(bestSubItr).instanceCenterIdx) - 1;
               remainingInstanceIdx = (IA(IA>= instanceOffset & IA <= instanceEndOffset) - instanceOffset) + 1;
               numberOfInstances = numel(remainingInstanceIdx);
               labelIds = num2cell(repmat(int32(remainingSubItr), numberOfInstances,1));

               %% Create required fields such as centerIdx, edges, children and sign array to assign to instances.
               centerIdx = bestSubs(bestSubItr).instanceCenterIdx(remainingInstanceIdx,:);
               matchCosts = bestSubs(bestSubItr).instanceMatchCosts(remainingInstanceIdx,:);
               adaptiveThreshold = single(optimalThreshold * ((size(bestSubs(bestSubItr).edges,1)+1)*2-1)) + singlePrecision;
               % Get graph match costs.
               instanceMatchScores = fix(multiplier * ((adaptiveThreshold - matchCosts)/ adaptiveThreshold)) / multiplier;
               numberOfInstances = numel(centerIdx);
               
               % Get instance children
               instanceChildren = instanceChildrenDescriptors(actualInstanceOffset:(actualInstanceOffset + numberOfInstances-1),:);
               instanceChildren( :, all(~instanceChildren,1) ) = [];
               
               % Calculate activations for the next level.
               instancePrevActivations = prevActivations(instanceChildren);
               if numberOfInstances>1
                   instanceMeanActivations = mean(instancePrevActivations,2);
               else
                   instanceMeanActivations = mean(instancePrevActivations);
               end
               instanceActivations = num2cell(instanceMatchScores .* instanceMeanActivations);
               instanceChildren = mat2cell(instanceChildren, ones(numberOfInstances,1), size(instanceChildren,2));
               instanceSigns = num2cell(allSigns(centerIdx));

               % Assign fields to graphs.
               [graphLevel(actualInstanceOffset:(actualInstanceOffset + numberOfInstances-1)).labelId] = deal(labelIds{:});
               [graphLevel(actualInstanceOffset:(actualInstanceOffset + numberOfInstances-1)).activation] = deal(instanceActivations{:});
               [graphLevel(actualInstanceOffset:(actualInstanceOffset + numberOfInstances-1)).children] = deal(instanceChildren{:});
               [graphLevel(actualInstanceOffset:(actualInstanceOffset + numberOfInstances-1)).sign] = deal(instanceSigns{:});
               instanceOffset = instanceEndOffset + 1;
               actualInstanceOffset = actualInstanceOffset + numberOfInstances;
               remainingSubItr = remainingSubItr + 1;
           end
        end
    else
        vocabLevel = [];
        graphLevel = [];
    end
   
    
   clear allEdges allEdgeNodePairs allSigns bestSubs remainingInstanceLabels instanceChildrenDescriptors;
   nextVocabLevel = vocabLevel;
   nextGraphLevel = graphLevel;
   clear vocabLevel graphLevel;
end

%> Name: getSingleNodeSubs
%>
%> Description: getSingleNodeSubs is used to obtain single-node subs from
%> labels of all instances in allLabels. The result is a number of
%> subs representing compositions each having their instances.
%> 
%> @param allLabels Labels for every graph node.
%> @param allSigns Signs for every graph node.
%> @param categoryArrIdx Categories for every graph node.
%> @param validationIdx Binary array showing if each node is part 
%> of the validation data or not.
%>
%> @retval singleNodeSubs Substructure list of single-node subs.
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 24.02.2014
%> Ver 1.1 on 01.09.2014 Removal of global parameters.
function singleNodeSubs = getSingleNodeSubs(allLabels, allSigns, nodeDistanceMatrix, categoryArrIdx, validationIdx, threshold)
    numberOfSubs = max(allLabels);
    singleNodeSubs(numberOfSubs) = Substructure();
    validSubs = ones(numberOfSubs,1)>0;
    
    %% For each center node label type, we create a substructure.
    for subItr = 1:numberOfSubs
        distances = nodeDistanceMatrix(allLabels, subItr);
        subCenterIdx = distances < threshold;
        instances = int32(find(subCenterIdx));
        numberOfInstances = numel(instances);
        
        %Assign center id.
        singleNodeSubs(subItr).centerId = subItr;
        
        % Give maximum score so that it is at the top of the queue.
        singleNodeSubs(subItr).mdlScore = numberOfInstances;
        if numberOfInstances>0
            % Fill in instance information. 
            categoryIdx = categoryArrIdx(subCenterIdx);
            if size(categoryIdx,1) == 1
                categoryIdx = categoryIdx';
            end
            instanceValidationIdx = validationIdx(subCenterIdx);
            instanceSigns = allSigns(subCenterIdx,1);
            instanceDistances = distances(subCenterIdx);
            
            % We check if this sub has exact-matching instances in the
            % training set. If not, it is not considered for further expansion.
            if nnz(instanceDistances == 0) == 0
                 singleNodeSubs(subItr).mdlScore = 0;
                 continue;
            end
            
            % Assign fields of the sub.
            singleNodeSubs(subItr).instanceCenterIdx = instances;
            singleNodeSubs(subItr).instanceChildren = instances;
            singleNodeSubs(subItr).instanceCategories = categoryIdx;
            singleNodeSubs(subItr).instanceMatchCosts = instanceDistances;
            singleNodeSubs(subItr).instanceSigns = instanceSigns;
            singleNodeSubs(subItr).instanceValidationIdx = instanceValidationIdx;
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
function [extendedSubs] = extendSub(sub, allEdges, nodeDistanceMatrix, edgeDistanceMatrix, singlePrecision, threshold)
    
    centerIdx = sub.instanceCenterIdx;
    subAllEdges = {allEdges(centerIdx).adjInfo}';
    % Get unused edges from sub's instances.
    allUsedEdgeIdx = sub.instanceEdges;
    if ~isempty(allUsedEdgeIdx)
        allUsedEdgeIdx = mat2cell(allUsedEdgeIdx, ones(size(allUsedEdgeIdx,1),1), size(allUsedEdgeIdx,2));
        allUnusedEdgeIdx = cellfun(@(x,y) setdiff(1:size(x,1),y), subAllEdges, allUsedEdgeIdx, 'UniformOutput', false);
        allUnusedEdges = cellfun(@(x,y) x(y,:), subAllEdges, allUnusedEdgeIdx, 'UniformOutput', false);
        allUnusedEdgeIdx = [allUnusedEdgeIdx{:}]';
    else
        allUnusedEdges = subAllEdges;
        allUnusedEdgeIdx = [];
    end
    
    % Record which edge belongs to which instance. 
    allEdgeInstanceIds = zeros(sum(cellfun(@(x) size(x,1), allUnusedEdges)),1);
    allEdgePrevCosts = zeros(size(allEdgeInstanceIds), 'single');
    itrOffset = 1;
    unusedEdgeCounts = cellfun(@(x) size(x,1), allUnusedEdges);
    for itr = 1:numel(allUnusedEdges)
        beginOffset = itrOffset;
        endOffset = (beginOffset+(unusedEdgeCounts(itr)-1));
        allEdgeInstanceIds(beginOffset:endOffset) = itr;
        allEdgePrevCosts(beginOffset:endOffset) = sub.instanceMatchCosts(itr);
        itrOffset = itrOffset + unusedEdgeCounts(itr);
    end         
    allUnusedEdges = cat(1, allUnusedEdges{:});
    
    % If no edges remain, exit.
    if isempty(allUnusedEdges)
        extendedSubs = [];
        return; 
    end
    
    % Eliminate the edges which exist only in validation data. We do not
    % enumerate any edges which do not exist in training data.
    enumeratedEdges = allUnusedEdges(allEdgePrevCosts < singlePrecision, 3:4);
    
    % Get unique rows of [edgeLabel, secondVertexLabel]
    uniqueEdgeTypes = unique(enumeratedEdges, 'rows');
    
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
    
    % Save local indices for all edges, to be used later.
    allLocalEdgeIdx = cellfun(@(x) 1:size(x,1), subAllEdges, 'UniformOutput', false);
    allLocalEdgeIdx = [allLocalEdgeIdx{:}]';
    
    % Allocate space for new subs and fill them in.
    extendedSubs(numberOfEdgeTypes) = Substructure;
    for edgeTypeItr = 1:numberOfEdgeTypes
        % Assign sub-definition type and other info
        newSub = Substructure;
        newSub.edges = sub.edges;
        newSub.centerId = sub.centerId;
        newSub.mdlScore = 0;
        newSub.edges = [newSub.edges; uniqueEdgeTypes(edgeTypeItr,:)];
        
        %% Find instances of this new sub, and mark the new edges as 'used' in each of its subs.
        % Extend the subs by taking the threshold into account. Unless the
        % instance's collective matching cost surpasses the threshold, it
        % is a valid instance.
        edgesToExtendCosts = allEdgePrevCosts + ...
            edgeDistanceMatrix(allUnusedEdges(:,3), uniqueEdgeTypes(edgeTypeItr,1)) + ...
            nodeDistanceMatrix(allUnusedEdges(:,4), uniqueEdgeTypes(edgeTypeItr,2));
        edgesToExtendIdx = edgesToExtendCosts < threshold;
        
        % Save instance ids.
        edgeInstanceIds = allEdgeInstanceIds(edgesToExtendIdx);
        allChildren = [sub.instanceChildren(edgeInstanceIds,:), ...
             allUnusedEdges(edgesToExtendIdx,2)];
        allChildren = sort(allChildren, 2);
        
        %% Note: allChildren may have duplicate entries, since an instance can
        % have multiple parse trees. We handle these cases by only keeping
        % unique instances. In addition, for each instance, the minimum
        % cost of matching is kept here.
        curInstanceMatchCosts = edgesToExtendCosts(edgesToExtendIdx);
        [minMatchCosts, sortIdx] = sort(curInstanceMatchCosts, 'ascend');
        sortedAllChildren = allChildren(sortIdx, :);
        sortedCenterIdx = sub.instanceCenterIdx(edgeInstanceIds(sortIdx));

        % Eliminate instances which have the same node set, and the same
        % center node. 
        [~, validIdx, ~] = unique([sortedCenterIdx, sortedAllChildren], 'rows', 'stable');
        
        % Keep ordering but remove duplicates.
        sortIdx = sortIdx(validIdx);
        
        % Get minimum matching costs and children.
        minMatchCosts = minMatchCosts(validIdx, :);
        sortedAllChildren = sortedAllChildren(validIdx, :);
        
        % Finally, order children by rows.
        [newSub.instanceChildren, idx] = sortrows(sortedAllChildren);
        sortIdx = sortIdx(idx);
        edgeInstanceIds = edgeInstanceIds(sortIdx, :);
        
        %% Assign all relevant instance-related fields of the sub.
        newSub.instanceCenterIdx = sub.instanceCenterIdx(edgeInstanceIds);
        newSub.instanceSigns = sub.instanceSigns(edgeInstanceIds);
        newSub.instanceCategories = sub.instanceCategories(edgeInstanceIds);
        newSub.instanceMatchCosts = minMatchCosts(idx,:);
        newSub.instanceValidationIdx = sub.instanceValidationIdx(edgeInstanceIds);
        
        % Add the edge to the definition.
        existingEdges = sub.instanceEdges;
        combinedEdges = zeros(numel(edgeInstanceIds), size(existingEdges,2)+1, 'uint8');
        if ~isempty(existingEdges)
            existingEdges = existingEdges(edgeInstanceIds,:);
            combinedEdges(:, 1:size(existingEdges,2)) = existingEdges;
        end
        
        % Adding the local edges to the instance edge lists.
        edgesToExtendLinIdx = find(edgesToExtendIdx);
        if ~isempty(sub.edges)
            relevantLocalEdgeIdx = allUnusedEdgeIdx(edgesToExtendLinIdx(sortIdx));
        else
            relevantLocalEdgeIdx = allLocalEdgeIdx(edgesToExtendLinIdx(sortIdx));
        end
        combinedEdges(:,end) = relevantLocalEdgeIdx;
        newSub.instanceEdges = combinedEdges;
        
        % All instances assigned, good to go.
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
function [subs] = evaluateSubs(subs, evalMetric, allEdges, allEdgeNodePairs, allSigns, graphSize, overlap, mdlNodeWeight, mdlEdgeWeight, isMDLExact, isSupervised)
    numberOfSubs = numel(subs);
    parfor subItr = 1:numberOfSubs
        % Find the weight of this node, by taking the max of the category distribution. 
        if isSupervised
            categoryArr = double(subs(subItr).instanceCategories);
            weight = nnz(categoryArr == mode(categoryArr)) / numel(categoryArr);
        else
            weight = 1;
        end
        
        % We compress the object graph using the children, and the
        % edges they are involved. 
        [subScore, sub] = getSubScore(subs(subItr), allEdges, allEdgeNodePairs, evalMetric, ...
           allSigns, mdlNodeWeight, mdlEdgeWeight, ....
            overlap, isMDLExact);
        subScore = subScore * weight;
        subs(subItr) = sub;

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
    clear sortedQueue addedQueue;
end

%> Name: removeDuplicateSubs
%>
%> Description: Given the child subs in childSubArr, this function removes
%> duplicate substructures from childSubArr, and combines the instances of
%> matching subs. If an instance can be parsed in different ways (i.e.
%> matches multiple duplicate subs), the minimum matching cost of matching is
%> saved. 
%> 
%> @param childSubArr A list of children substructures.
%>
%> @retval childSubArr The list of unique substructures.
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 01.03.2015 (Converted to a standalone function, added instance augmentation)
function [childSubArr] = removeDuplicateSubs(childSubArr)
    if size(childSubArr(1).edges,1) > 1
        % Eliminate duplicate subs in final array.
        childSubEdges = cell(1, numel(childSubArr));
        for subItr = 1:numel(childSubArr)
            childSubEdges(subItr) = {sortrows(childSubArr(subItr).edges)};
        end
        subCenters = {childSubArr.centerId};
        vocabDescriptors = cellfun(@(x,y) num2str([x; y(:)]'), subCenters, childSubEdges, 'UniformOutput', false);
        [~, validChildrenIdx, ~] = unique(vocabDescriptors, 'stable');
        childSubArr = childSubArr(validChildrenIdx);
    end
end