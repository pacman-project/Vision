% this function provides a unified interface to all functions reading
% diffrent dataSets!

function [list_depth, list_mask, list_RGB, lenF] = extractFileListGeneral(depthPath, is_subset, subsetPercent, dataSetNumber)


    if dataSetNumber == 1 || dataSetNumber == 3
        
        [list_depth, lenF] = extractFileList(depthPath, is_subset, subsetPercent);
        list_mask = [];
        list_RGB = [];
        
    elseif dataSetNumber == 2
        
        [list_depth, list_mask, list_RGB, lenF] = extractFileListWashington(depthPath, is_subset, subsetPercent);
        
    elseif dataSetNumber == 4  % for Mirela Dataset
        
        
        is_RGB = false;
        is_Mask = false;
        [list_depth, list_mask, list_RGB, lenF] = extractFileListMirela(depthPath, is_subset, subsetPercent, is_RGB, is_Mask);
        
    elseif dataSetNumber == 5  % DataSet of Meshes of Aim&Shape
        
        [list_depth, lenF] = extractFileList_meshes(depthPath, is_subset, subsetPercent);
        list_mask = [];
        list_RGB = [];
        
    end

end

