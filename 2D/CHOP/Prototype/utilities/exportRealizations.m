%> Name: exportRealizations
%>
%> Description: Export all realizations in all object graph layers to a
%> single matrix 'exportArr'. It is of size N x 5, where N = number of
%> realizations. 
%>
%> @param mainGraph Object graphs.
%>
%> @retval exportArr Matrix including information for all realizations,
%> of the form: 
%> labelId1 coordX1 coordY1 levelId1 imageId1;
%> labelId2 coordX2 coordY2 levelId2 imageId2;
%> ...
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 19.12.2013
function [exportArr, activationArr, pooledPositions] = exportRealizations(graphLevel, levelItr)
    %% Process mainGraph to export realizations in the desired format for inte2D/3D integration.
    % n x 5 array. Each row represents a realization. Format is:
    % labelId coordX coordY levelId imageId.
    labelIds = [graphLevel.realLabelId]';
    activationArr = [graphLevel.activation]';
    positions = cat(1, graphLevel.precisePosition);
    pooledPositions = cat(1, graphLevel.position);
    imageIds = [graphLevel.imageId]';
    levelIds = ones(numel(labelIds),1, 'int32') * levelItr;
    exportArr = [labelIds, positions, levelIds, imageIds];
end