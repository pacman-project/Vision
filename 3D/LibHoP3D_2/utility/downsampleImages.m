% this function downsamples depth images and mask images according to the size of MARKS image 

function [ ] = downsampleImages(list_depths, depthPath, list_El, list_mask,  dataSetNumber, outFolder)

    disp('Downsampling of images and masks...');
    
    lenF = length(list_depths);
    lenDPW = length(depthPath);
    if isempty(list_mask)
        list_mask = zeros(1,lenF);
    end
    
    for i = 1:lenF
        I = imread(list_depths{i});
        I = I(:,:,1);
        
        marks = imread(list_El{i});
        [r,c] = size(marks);
      
        if dataSetNumber == 2
            
            mask = imread(list_mask{i});
            mask = imresize(mask, [r, c]);

        end
        
        I = imresize(I, [r, c]);

        
        % then we store the output image (or images to the folder outFolder)
        
        curStr = list_depths{i};
        fileName = curStr(lenDPW+1:end);
        outFile = [outFolder, fileName];

        ll = strfind(outFile, '/');
        ll = ll(end); % last position
        folderName = outFile(1:ll);
        b = exist(folderName,'dir');

        if b == 0
            mkdir(folderName);
        end
        
        if dataSetNumber == 2 % we also need to store a downsampled mask
            
            curStr = list_mask{i};

            fileName = curStr(lenDPW+1:end);
            outFile = [outFolder, fileName];

%             ll = strfind(outFile, '/');
%             ll = ll(end); % last position
%             folderName = outFile(1:ll);
%             b = exist(folderName,'dir');
% 
%             if b == 0
%                 mkdir(folderName);
%             end
        end
        
        if dataSetNumber == 2
            disp('Error');
        end

        imwrite(I, outFile, 'png');
        
        if mod(i,20) == 0
            i
        end
    end
end
