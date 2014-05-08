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
function [graphLevel] = inferSubs(vocabLevel, graphLevel, options)
    global redundantVocabLevel;
    % Combine vocabLevel with redundant compositions.
    validVocabLevel = [ones(numel(vocabLevel),1); zeros(numel(redundantVocabLevel),1)];
    realVocabIds = [1:numel(vocabLevel), [redundantVocabLevel.label]];
    vocabLevel = [vocabLevel, redundantVocabLevel];
    
    % Read data into helper data structures.
    edges = {graphLevel.adjInfo}';
    
    % If no edges are present, time to return.
    if isempty(edges)
        graphLevel = [];
        return;
    end
    
    %% Read edges and label ids of linked nodes into edges array.
    nonemptyEdgeIdx = cellfun(@(x) ~isempty(x), edges);
    labelIds = cat(1, graphLevel.labelId)';
    edges(nonemptyEdgeIdx) = cellfun(@(x) [x(:,1:3), labelIds(x(:,1:2))], ...
        edges(nonemptyEdgeIdx), 'UniformOutput', false);
    edgeCount = cellfun(@(x) size(x,1), edges);
    edges = cat(1, edges{nonemptyEdgeIdx});
    
    % If no edges are present, time to return.
    if isempty(edges)
        graphLevel = [];
        return;
    end
    
    % Record which edge belongs to which instance. 
    allEdgeInstanceIds = zeros(size(edges,1),1);
    itrOffset = 1;
    for itr = 1:numel(edgeCount)
        beginOffset = itrOffset;
        allEdgeInstanceIds(beginOffset:(beginOffset+(edgeCount(itr)-1))) = itr;
        itrOffset = itrOffset + edgeCount(itr);
    end 
    
    %% Get descriptors for edges in the object graph.
    if strcmp(options.property, 'mode')
        edgeDescriptors = edges(:,3);
    else
        offset1 = max(max(edges(:,4:5)));
        offset2 = offset1^2;
        edgeDescriptors = (edges(:,3)-1) * offset2 + (edges(:,4)-1)*offset1 + edges(:,5);
    end
    
    %% Get unique label ids of the previous level.
    uniqueLabelIds = unique(labelIds);
    
    %% Match subs from vocabLevel to their instance in graphLevel.
    vocabRealizations = cell(numel(vocabLevel),1);
    for vocabItr = 1:numel(vocabLevel)
       %% Quick checks to determine if this composition should be searched.
       % First, see if all nodes exist in the object graphs.
       if nnz(ismember(vocabLevel(vocabItr).children, uniqueLabelIds)) ~= numel(vocabLevel(vocabItr).children)
           continue;
       end
       
       %% Get descriptors for edges in the vocabulary node.
       vocabEdges = vocabLevel(vocabItr).adjInfo;
       vocabChildren = (vocabLevel(vocabItr).children)';
       if strcmp(options.property, 'mode')
            vocabDescriptors = vocabEdges(:,3);
       else
            vocabDescriptors = (vocabEdges(:,3)-1)*offset2 + (vocabChildren(vocabEdges(:,1))-1) * offset1 + vocabChildren(vocabEdges(:,2));
       end
       numberOfVocabEdges = numel(vocabDescriptors);
       
       %% INDEXING: Search for descriptors in receptive fields, and return those contain at least one.
       nonUnqValidInstances = allEdgeInstanceIds(edgeDescriptors == vocabDescriptors(1));
       if ~isempty(nonUnqValidInstances)
           validInstances = fastsortedunique(allEdgeInstanceIds(edgeDescriptors == vocabDescriptors(1)));
       else
           validInstances = [];
       end
       
       % Verify valid instances by checking against rest of edges.
       if numberOfVocabEdges > 1 && ~isempty(validInstances)
           for vocabDescItr = 2:numberOfVocabEdges
                validInstances = validInstances(ismembc(validInstances, ...
                    allEdgeInstanceIds(edgeDescriptors == vocabDescriptors(vocabDescItr))));
           end
       end    
       
       % If no valid instances exist, pass.
       if isempty(validInstances)
           continue;
       end
       
       %% Process each instance to get correct node sets corresponding to instances.
       discoveredInstances = cell(numel(validInstances),1);
       for instanceItr = 1:numel(validInstances)
            instanceEdgeIdx = find(allEdgeInstanceIds == validInstances(instanceItr));
            instanceEdgeDescriptors = edgeDescriptors( allEdgeInstanceIds == validInstances(instanceItr));
            
            % Collect instances of each edge from vocabulary node's description.
            instanceEdgeSets = cell(numberOfVocabEdges,1);
            for vocabEdgeItr = 1:numberOfVocabEdges
                idx = instanceEdgeDescriptors == vocabDescriptors(vocabEdgeItr);
                instanceEdgeSets(vocabEdgeItr) = {edges(instanceEdgeIdx(idx), 2)};
            end
            
            % Get all possible combinations of nodes from each set.
            discoveredInstances(instanceItr) = {allcomb(validInstances(instanceItr), instanceEdgeSets{:})};
       end
       vocabRealizations(vocabItr) = {cat(1, discoveredInstances{:})};
    end
    numberOfInstanceArr = cellfun(@(x) size(x,1), vocabRealizations);
    numberOfInstances = sum(numberOfInstanceArr);
    clear graphLevel;
    
    %% If no instances have been found, exit.
    if numberOfInstances<1
        graphLevel = [];
        return;
    end
    
    %% Generate graph level.
    graphLevel(numberOfInstances) = options.graphNode;
    
    % Collect label ids, children and signs to assign to instances.
    labelIds = zeros(numberOfInstances,1);
    instanceOffset = 1;
    for vocabItr = 1:numel(vocabLevel)
       if numberOfInstanceArr(vocabItr) > 0
            instanceEndOffset = instanceOffset + (numberOfInstanceArr(vocabItr)-1);
            if validVocabLevel(vocabItr)
                labelIds(instanceOffset:instanceEndOffset) = vocabItr;
            else
                labelIds(instanceOffset:instanceEndOffset) = realVocabIds(vocabItr);
            end
            instanceOffset = instanceEndOffset+1;
       end
    end
    labelIds = num2cell(labelIds);
    allChildren = cellfun(@(x) mat2cell(x, ones(size(x,1),1), size(x,2)), vocabRealizations, 'UniformOutput', false);
    allChildren = cat(1, allChildren{:});
    
    %% Assign labelId, children and sign of the children.
    [graphLevel.labelId] = deal(labelIds{:});
    [graphLevel.children] = deal(allChildren{:});
    [graphLevel.sign] = deal(1);
end

