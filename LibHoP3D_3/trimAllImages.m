
fileListPrecomputed = false;
is_subset = false;
dataSetNumber = 1;
subset_len = 400;

depthPathDefault = '';
depthPath = 'D:/Input Data/AimShape/4T/';  
ldp = length(depthPath);
isTrim = true;

outFolder = 'D:/Input Data/AimShape/4TT/';  %'D:/3D/Input Data/Images for categorization/1T/';

if dataSetNumber == 1 || dataSetNumber == 3
    [list_depth, lenF] = extractFileList(fileListPrecomputed, depthPath, depthPathDefault, is_subset, subset_len);
    list_mask = [];
elseif dataSetNumber == 2
    [list_depth, list_mask, ~, lenF] = extractFileListWashington(fileListPrecomputed, depthPath, depthPathDefault, is_subset, subsetPercent);
end

parfor i = 1 : lenF
    
    str = list_depth{i};
    I = imread(str);    

    [r,c] = size(I);
    mask = zeros(r,c);
    mask(I > 0) = 1;

    % check if the mask is empty or not
    maxM = max(max(mask));
    if maxM == 0 
        is_successfull = false;
    end
    if isTrim
        [I, mask] = trimImage(I, mask); % trim the image
        
        % save the trimmed image into the new folder
        fileName = str(ldp+2:end);
        outFileName = [outFolder, fileName];
        [r,c] = size(I);
        I3 = zeros(r,c,3);
        for k = 1:3
            I3(:,:,k) = I; % all channels are the same
        end
        I3 = uint16(I3);
        imwrite(I3, outFileName, 'png');
        
    end
    
    
end
    
    