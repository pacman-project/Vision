% this function reproduces all images at all required scales

function computeAllscales(scales, inputPath, lineAdders, dataSetNumber, inputDataType)

    % check if downsampled folders exist
    if inputDataType == 1 % depth image
    
        for  i = 1:length(scales)
            curFolder = [inputPath, '_', lineAdders{i}];

            if ~exist(curFolder, 'dir') 
                copyfile(inputPath, curFolder);
                UpsampleAllImages(dataSetNumber, curFolder, scales(i), false, 1.0);
            end
        end
    elseif inputDataType == 2 % meshes
        
        for  i = 1:length(scales)
            curFolder = [inputPath, '_', lineAdders{i}];

            if ~exist(curFolder, 'dir') 
                copyfile(inputPath, curFolder);
                UpsampleAllMeshes(dataSetNumber, curFolder, scales(i), false, 1.0);
            end
        end
        
    end
    
    

end

