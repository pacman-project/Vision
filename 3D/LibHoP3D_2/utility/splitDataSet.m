% this is a function to split a Washington data set

dataSetNumber = 2; % Washington

depthPath_W = '/home/vvk201/Wash-rgbd-dataset';
depthPathDefault = '';

outputFolder = '/home/vvk201/Wash-rgbd-dataset_0003';

fileListPrecomputed = false;
is_subset = true;
subsetPercent = 0.003;


if dataSetNumber == 1
    [list_depth, lenF] = extractFileList(fileListPrecomputed, depthPath, depthPathDefault, is_subset, subset_len);
    list_mask = [];
elseif dataSetNumber == 2
    [list_depth, list_mask, list_images, lenF] = extractFileListWashington(fileListPrecomputed, depthPath_W, depthPathDefault, is_subset, subsetPercent);
end

strFolderLen = length(depthPath_W);


parfor i = 1:lenF % we have to copy each file to the new location
    
    strD = list_depth{i};
    strI = list_images{i};
    strM = list_mask{i};
    
    newFileD = [outputFolder, strD(strFolderLen + 1:end)];
    newFileI = [outputFolder, strI(strFolderLen + 1:end)];
    newFileM = [outputFolder, strM(strFolderLen + 1:end)];
    
    ll = strfind(newFileD, '/');
    ll = ll(end);
    
    curFolderName = newFileD(1:ll);
    A = exist(curFolderName, 'dir');
    
    if A == 0
        mkdir(curFolderName);
    end
    
    copyfile(list_depth{i}, newFileD);
    copyfile(list_images{i}, newFileI);
    copyfile(list_mask{i}, newFileM);
end