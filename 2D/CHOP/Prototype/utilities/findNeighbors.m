%> Name: findNeighbors
%>
%> Description: This function finds 1- and 2-neighborhood adjacent nodes
%> for each first layer node. 
%>
%> @param firstLevelAdjNodes Adjacency information (immediate) for the
%> first layer. Format for every row: [centerId, adj1, adj2, ..., 0, 0, 0].
%>
%> @retval firstLevelAdjNodes1 1-neighbors. Same format as input.
%> @retval firstLevelAdjNodes2 2-neighbors. Same format as input.
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 05.01.2016
function [firstLevelAdjNodes1, firstLevelAdjNodes2] = findNeighbors(firstLevelAdjNodes)
    firstLevelAdjNodes1 = cell(size(firstLevelAdjNodes,1),1);
    firstLevelAdjNodes2 = cell(size(firstLevelAdjNodes,1),1);
    max1 = 0; max2 = 0;
    for nodeItr = 1:size(firstLevelAdjNodes,1)
        % Obtain immediate neighbors.
        curNeighbors = firstLevelAdjNodes(nodeItr,:);
        curNeighbors = sort(curNeighbors(curNeighbors>0));
        neighbors = curNeighbors(2:end);
        neighbors = neighbors(neighbors>0);
        
        % No neighbors? Move on.
        if nnz(neighbors) == 0
           continue; 
        end
        
        % Obtain indirect neighbors (1-neigh)
        neighbors1 = firstLevelAdjNodes(neighbors,2:end);
        neighbors1 = neighbors1(:)';
        
        % No 1-neighbors? Move on.
        if nnz(neighbors1) == 0
           continue; 
        end
        
        % Save 1-neighbors.
        neighbors1 = fastsortedunique(sort(neighbors1(neighbors1>0)));
        neighbors1 = neighbors1(~ismembc(neighbors1, curNeighbors));
        firstLevelAdjNodes1{nodeItr} = neighbors1;
        max1 = max(max1, numel(neighbors1));
        
        % Obtain indirect neighbors (2-neigh)
        neighbors2 = firstLevelAdjNodes(neighbors1,2:end);
        neighbors2 = neighbors2(:)';
        
        % No 1-neighbors? Move on.
        if nnz(neighbors2) == 0
           continue; 
        end
        
        curNeighbors = sort(cat(2, neighbors1, curNeighbors));
        
        % Save 2-neighbors.
        neighbors2 = fastsortedunique(sort(neighbors2(neighbors2>0)));
        neighbors2 = neighbors2(~ismembc(neighbors2, curNeighbors));
        firstLevelAdjNodes2{nodeItr} = neighbors2;
        max2 = max(max2, numel(neighbors2));
    end
    % Process the neighborhood arrays so they have the same format as the
    % first one.
    % 1-neighbors
    dummyArrs = cell(max1,1);
    for itr = 1:max1
        dummyArrs(itr) = {zeros(1, max1 - itr, 'int32')};
    end
    nonzeroIdx = cellfun(@(x) ~isempty(x), firstLevelAdjNodes1);
    firstLevelAdjNodes1(nonzeroIdx) = cellfun(@(x) [x, dummyArrs{numel(x)}], firstLevelAdjNodes1(nonzeroIdx), 'UniformOutput', false);
    firstLevelAdjNodes1(~nonzeroIdx) = {zeros(1, max1, 'int32')};
    firstLevelAdjNodes1 = cat(1, firstLevelAdjNodes1{:});
    
    % 2-neighbors
    dummyArrs = cell(max2,1);
    for itr = 1:max2
        dummyArrs(itr) = {zeros(1, max2 - itr, 'int32')};
    end
    nonzeroIdx = cellfun(@(x) ~isempty(x), firstLevelAdjNodes2);
    firstLevelAdjNodes2(nonzeroIdx) = cellfun(@(x) [x, dummyArrs{numel(x)}], firstLevelAdjNodes2(nonzeroIdx), 'UniformOutput', false);
    firstLevelAdjNodes2(~nonzeroIdx) = {zeros(1, max2, 'int32')};
    firstLevelAdjNodes2 = cat(1, firstLevelAdjNodes2{:});
end