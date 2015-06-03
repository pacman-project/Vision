% this function mearges all layers to one multichannel image

function [] = meargeAllLayers()

    root = '/home/vvk201/LibHoP3D/'; % 'D:/3D/LibHoP3D/'; 
    addPaths(root);

    nClusters = 7;
    n2Clusters = nClusters^2;
    dataSetNumber = 2;
    abstractionLevel = 3;
    is_subset = false;
    
    dsN = num2str(dataSetNumber);
    nCl = num2str(nClusters);
    aL = num2str(abstractionLevel);

    
    parts3Layer = [root, 'statistics/partsSelectionResults_3_', dsN, '_', nCl, '_a', aL, '.mat'];
    parts4Layer = [root, 'statistics/partsSelectionResults_4_', dsN, '_', nCl, '_a', aL, '.mat'];
    parts5Layer = [root, 'statistics/partsSelectionResults_5_', dsN, '_', nCl, '_a', aL, '.mat'];
    parts6Layer = [root, 'statistics/partsSelectionResults_6_', dsN, '_', nCl, '_a', aL, '.mat'];
    
    % define paths here
    
    pathsLayers{1} = '/home/vvk201/Wash-rgbd-dataset_layer2';
    pathsLayers{2} = '/home/vvk201/Wash-rgbd-dataset_layer3';
    pathsLayers{3} = '/home/vvk201/Wash-rgbd-dataset_layer4';
    
    lenDPW = length(pathsLayers{1});
    
    outPath = '/home/vvk201/Wash-rgbd-dataset_layers234';
    
%     load(parts3Layer);  % to define n3Clusters and n4Clusters
%     load(parts4Layer);

    n3Clusters = 150;
    n4Clusters = 320;

    if dataSetNumber == 2  % extract images with previous layer
        [list_El2, ~, ~, lenF2] = extractFileListWashingtonForClassification(pathsLayers{1}, is_subset, 1.0);
        [list_El3, ~, ~, lenF3] = extractFileListWashingtonForClassification(pathsLayers{2}, is_subset, 1.0);
        [list_El4, ~, ~, lenF4] = extractFileListWashingtonForClassification(pathsLayers{3}, is_subset, 1.0);
    else
        disp('ERROR');
        return;
    end
    
    if lenF2 ~= lenF3 || lenF2 ~= lenF4
        disp('Error');
        return;
    end
    
    parfor i = 1:lenF2
        
        marks2 = imread(list_El2{i});
        marks3 = imread(list_El3{i});
        marks4 = imread(list_El4{i});
        
        [r,c] = size(marks2);
        I = zeros(r,c,3);
        I = uint16(I);
        
    
        I(:,:,1) = uint16(marks2);
        
        shift = n2Clusters + 1;
        marks3 = marks3 + shift;
        marks3(marks3 == shift) = 0;
        I(:,:,2) = marks3;
        
        shift = n2Clusters + n3Clusters + 1;
        marks4 = marks4 + shift;
        marks4(marks4 == shift) = 0;
        I(:,:,3) = marks4;
        
  
        % save the image I
        
        curStr = list_El2{i};

        fileName = curStr(lenDPW+1:end);
        outFile = [outPath, fileName];

        ll = strfind(outFile, '/');
        ll = ll(end); % last position
        folderName = outFile(1:ll);
        b = exist(folderName,'dir');

        if b == 0
            mkdir(folderName);
        end

        I = uint16(I);
        imwrite(I, outFile, 'png');
    end
        
    
    
end