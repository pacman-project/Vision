function UpsampleAllImages(dataSetNumber, depthPath, dowsample_rate)

    fileListPrecomputed = false;
    is_subset = false; % whether we shall use all files for learning

    % define the subset length
    if dataSetNumber == 1 || dataSetNumber == 3
        subset_len = 800; % how much shall we use for training
        subsetPercent = 1.0; % not used
    elseif dataSetNumber == 2
        subset_len = 1;
        subsetPercent = 0.2; % what percent from each folder to use
    end
    
    
    if dataSetNumber == 1 || dataSetNumber == 3
        [list_depth, lenF] = extractFileList(depthPath{1}, is_subset, subset_len);
        list_mask = [];
    elseif dataSetNumber == 2
        [list_depth, list_mask, ~, lenF] = extractFileListWashington(depthPath{1}, is_subset, subsetPercent);
    end
    
    upsampleImages(list_depth, list_mask, lenF, true, dowsample_rate); % to be done only once!


end

