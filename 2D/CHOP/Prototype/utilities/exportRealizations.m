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
function [exportArr] = exportRealizations(mainGraph)
    %% Process mainGraph to export realizations in the desired format for inte2D/3D integration.
    % n x 5 array. Each row represents a realization. Format is:
    % labelId coordX coordY levelId imageId.
    allLevels = [mainGraph{:}];
    labelIds = [allLevels.labelId]';
    positions = cat(1, allLevels.position);
    imageIds = [allLevels.imageId]';
    levelIds = zeros(numel(labelIds),1);
    offset = 1;
    for levelItr = 1:numel(mainGraph)
        realCount = numel(mainGraph{levelItr});
        levelIds(offset:(offset + (realCount - 1))) = levelItr;
        offset = offset + realCount;
    end
    exportArr = [labelIds, positions, levelIds, imageIds];
end