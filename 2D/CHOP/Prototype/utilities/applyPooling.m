%> Name: applyPooling
%>
%> Description: This function applies max pooling to the nodes represented
%> by graphLevel. 
%>
%> @param graphLevel The current graph level.
%> @param poolDim Pooling dimension.
%> @param poolFlag Flag for pooling-based elimination.
%>
%> @retval updatedGraphLevel Remaining graph level nodes.
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 16.09.2015
function [ updatedGraphLevel ] = applyPooling( graphLevel, poolFlag )
     labelIds = cat(1, graphLevel.labelId);
     imageIds = cat(1, graphLevel.imageId);
     coords = cat(1, graphLevel.position);
     activations = cat(1, graphLevel.activation);
     
     % First, put the children into their allocated spaces.
%     children = {graphLevel.children};
%     maxSize = max(cellfun(@(x) size(x,2), children));
%      vocabRealizationsChildren = cellfun(@(x) ...
%         cat(2, x, zeros(size(x,1), maxSize - size(x,2), 'int32')), ...
%         children, 'UniformOutput', false);
%     vocabRealizationsChildren = cat(1, vocabRealizationsChildren{:});
     
     % Get unique nodes for each label, image, coord triplet.
     combinedArr = double([imageIds, labelIds, coords]);
%     combinedArr = double([imageIds, coords]);
     
%      % First, we order the nodes by labelIds and coords.
%      arrayToSort = [combinedArr, double(vocabRealizationsChildren)];
%      [~, sortIdx] = sortrows(arrayToSort);
%      graphLevel = graphLevel(sortIdx);
%      combinedArr = combinedArr(sortIdx, :);
%      activations = activations(sortIdx,:);
     
     % Sort combinedArr so that it is sorted by decreasing activations.
     [~, idx] = sort(activations, 'descend');
     combinedArr = combinedArr(idx,:);
     
     % Perform pooling if needed.
     if poolFlag
          [~, IA, ~] = unique(combinedArr, 'rows', 'stable');
     else
          IA = (1:size(combinedArr,1))';
     end
     
     % Save real indices and activations.
     idx = idx(IA);
     idx = sort(idx);
     activations = activations(idx);
     
     %% Create final graphLevel.
     updatedGraphLevel = graphLevel(idx);
     activations = num2cell(activations);
     [updatedGraphLevel.activation] = deal(activations{:});
     
    % Rearrange graph level so it is sorted by image id.
    arrayToSort = [[updatedGraphLevel.imageId]', [updatedGraphLevel.labelId]', cat(1, updatedGraphLevel.position)];
    [~, sortedIdx] = sortrows(arrayToSort);
    updatedGraphLevel = updatedGraphLevel(sortedIdx);
    clearvars -except updatedGraphLevel
end

