%> Name: applyPooling
%>
%> Description: This function applies max pooling to the nodes represented
%> by graphLevel. 
%>
%> @param graphLevel The current graph level.
%> @param poolDim Pooling dimension.
%>
%> @retval updatedGraphLevel Remaining graph level nodes.
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 16.09.2015
function [ updatedGraphLevel ] = applyPooling( graphLevel, poolDim )
     labelIds = cat(1, graphLevel.labelId);
     imageIds = cat(1, graphLevel.imageId);
     coords = cat(1, graphLevel.position);
     activations = cat(1, graphLevel.activation);
     
     % Get unique nodes for each label, image, coord triplet.
     combinedArr = [labelIds, imageIds, coords];
     
     % Sort combinedArr so that it is sorted by decreasing activations.
     [~, idx] = sort(activations, 'descend');
     combinedArr = combinedArr(idx,:);
     
     %% Downsample the coordinates (pooling), and then perform max operation.
     combinedArr(:,3:4) = floor(combinedArr(:,3:4) / poolDim);
     [~, IA, ~] = unique(combinedArr, 'rows', 'stable');
     
     % Save real indices and activations.
     idx = idx(IA);
     idx = sort(idx);
     activations = activations(idx);
     
     %% Create final graphLevel.
     updatedGraphLevel = graphLevel(idx);
     updatedPositions = coords(idx,:);
     updatedPositions = int32(floor((double(updatedPositions) - 1)/poolDim) + 1);
     updatedPositions = mat2cell(updatedPositions, ones(numel(updatedGraphLevel),1), 2);
     activations = num2cell(activations);
     [updatedGraphLevel.activation] = deal(activations{:});
     [updatedGraphLevel.position] = deal(updatedPositions{:});
end

