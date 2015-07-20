function UpsampleAllMeshes(dataSetNumber, depthPath, dowsample_rate, is_subset, subsetPercent)
    
    [list_input, ~, ~, lenF] = extractFileListGeneral(depthPath, is_subset, subsetPercent, dataSetNumber);
    
    upsampleMeshes(list_input, lenF, dowsample_rate); % downsamples images in place!
end

