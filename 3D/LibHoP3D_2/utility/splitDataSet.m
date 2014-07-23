% this is a function to split a Washington data set

select_depths = true;
select_elements = false;

dataSetNumber = 2;
depthPathDefault = '';

if select_depths    
    list_mask = [];
    list_images = [];
%     depthPath = 'D:\Input Data\AimShape\4T';
%     outputFolder = 'D:\Input Data\AimShape\4T_600';

    depthPath = 'D:\Input Data\Washington\Washington3Categories';
    outputFolder = 'D:\Input Data\Washington\Washington3Categories_008';
   
end

if select_elements
    elPath = '/home/vvk201/Input Data/Washington/Wash-rgbd-dataset_layer4';
    outputFolderE = '/home/vvk201/Input Data/Washington/Wash-rgbd-dataset_layer4_008';
end

fileListPrecomputed = false;
is_subset = true;
subset_len = 600;
subsetPercent = 0.08;

strFolderLen = length(depthPath);

if select_elements
    strFolderLenE = length(elPath);
end

if dataSetNumber == 1 || dataSetNumber == 3
    [list_depth, lenF] = extractFileList(fileListPrecomputed, depthPath, depthPathDefault, is_subset, subset_len);
    list_images = zeros(1, lenF);
    list_mask = zeros(1, lenF);
elseif dataSetNumber == 2
    [list_depth, list_mask, list_images, lenF] = extractFileListWashington(fileListPrecomputed, depthPath, depthPathDefault,...
        is_subset, subsetPercent);
end

if select_elements
    % select the same files from folder elPath

    list_els = list_depth;
    for i = 1:lenF
        str = list_els{i};
        str = [elPath ,str(strFolderLen + 1:end)];
        list_els{i} = str;
        a = 2;
    end
end


for i = 1:lenF % we have to copy each file to the new location
    
    strD = list_depth{i};
    newFileD = [outputFolder, strD(strFolderLen + 1:end)];
    

    if dataSetNumber == 2
        strI = list_images{i};
        strM = list_mask{i};
        newFileI = [outputFolder, strI(strFolderLen + 1:end)];
        newFileM = [outputFolder, strM(strFolderLen + 1:end)];
    end
    if select_elements
        strE = list_els{i};
        newFileE = [outputFolderE, strD(strFolderLen + 1:end)];
    end
    
    ll = strfind(newFileD, '/');
    ll = ll(end);
    
    curFolderName = newFileD(1:ll);
    A = exist(curFolderName, 'dir');
    
    if A == 0
        mkdir(curFolderName);
    end
    
    copyfile(list_depth{i}, newFileD);
    
    if dataSetNumber == 2
        copyfile(list_images{i}, newFileI);
        copyfile(list_mask{i}, newFileM);
    end
    
    
    if select_elements
        ll = strfind(newFileE, '/');
        ll = ll(end);

        curFolderNameE = newFileE(1:ll);
        A = exist(curFolderNameE, 'dir');

        if A == 0
            mkdir(curFolderNameE);
        end

        copyfile(list_els{i}, newFileE);
    end
    
end