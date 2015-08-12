function [ nodes ] = projectNode( nodes, vocabulary, inhibitionRadius )
    levelItr = nodes(1,4);
    nodes = single(nodes);
    inhibitionRadius = inhibitionRadius - 1;
    inhibitionRadiusSq = inhibitionRadius^2;
    
    %% First, we recursively backproject the nodes. 
    while levelItr > 1.001
        vocabLevel = vocabulary{levelItr};
        newNodes = cell(size(nodes,1),1);
        for nodeItr = 1:size(nodes,1)
            vocabNode = vocabLevel(nodes(nodeItr,1));
            newNodeSet = zeros(numel(vocabNode.children), 4, 'single');
            newNodeSet(:, 1) = single((vocabNode.children)');
            newNodeSet(:, 4) = levelItr-1;
            newNodeSet(:, 2:3) = (vocabNode.childrenPosMean) - repmat(nodes(nodeItr, 2:3), size(newNodeSet,1), 1);
            newNodes{nodeItr} = newNodeSet;
        end
        nodes = cat(1, newNodes{:});
        levelItr = levelItr - 1;
    end
    nodes = int32(round(nodes(:,1:3)));
    
    %% Then, we perform a simple inhibition process to remove overlapping level 1 instances.
    validNodes = ones(size(nodes,1),1) > 0;
    numberOfNodes = size(nodes,1);
    for nodeItr = 1:(numberOfNodes-1)
        if ~validNodes(nodeItr) 
           continue; 
        end
        remainingNodes = nodes((nodeItr+1):end, :);
        validRemainingNodes = sum((remainingNodes(:,2:3) - repmat(nodes(nodeItr,2:3), ...
            numberOfNodes - nodeItr, 1)).^2, 2) >= inhibitionRadiusSq;
        validNodes((nodeItr+1):end) = validNodes((nodeItr+1):end) & validRemainingNodes;
    end
    
    nodes = nodes(validNodes, :);
end

