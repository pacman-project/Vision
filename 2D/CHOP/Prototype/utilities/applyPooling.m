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
function [ updatedGraphLevel ] = applyPooling( graphLevel, poolFlag, labeledPooling )
     labelIds = cat(1, graphLevel.labelId);
     imageIds = cat(1, graphLevel.imageId);
     coords = cat(1, graphLevel.position);
     activations = cat(1, graphLevel.activation);
     
     % Get unique nodes for each label, image, coord triplet.
     if labeledPooling
          combinedArr = double([imageIds, labelIds, coords]);
     else
          combinedArr = double([imageIds, coords]);
     end
     
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

