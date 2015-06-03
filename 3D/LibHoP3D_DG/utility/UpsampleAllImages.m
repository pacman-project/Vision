function UpsampleAllImages(dataSetNumber, depthPath, dowsample_rate, is_subset, subsetPercent)
    
    [list_depth, list_mask, ~, lenF] = extractFileListGeneral(depthPath, is_subset, subsetPercent, dataSetNumber);
    checkImages2(list_depth, lenF);
    
    upsampleImages(list_depth, list_mask, lenF, true, dowsample_rate); % downsamples images in place!
end

