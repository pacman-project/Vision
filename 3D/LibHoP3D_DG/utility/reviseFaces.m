% this fucntion goes through all images and through all faces and checks a
% minimal and average likelihood in the neighborhood.  If it is above a
% certain value - we do not perform inference for the bigger scales at this
% point

function crossScaleStructure = reviseFaces(crossScaleStructure, list_input, lenFiles, receptiveFieldNext, strE, lenSc, curScale)

    if isempty(crossScaleStructure)
        for i = 1:lenFiles
            crossScaleStructure{i} = {};
        end
    end
    
    threshMin = 0.4;
    threshAvg = 0.60;
    kkMax = 2;
    
    for kk = 1:lenFiles
        
        fullFileName = list_input{kk};        
%         fileName = 'D:\Input Data\Meshes\Aim@Shape_Selected_2.00\D00002.obj';

        [folderName, fileName] = getFolderName(fullFileName);
        fileM =  [strE, '/', fileName(1:end-4), '.mat']; % file with parts
        fileAP = [strE, '/', fileName(1:end-4), 'AP.mat']; % file with additional points

        [~, F, ~] = meshRead(fullFileName);
        load(fileM); %  'Vout', 'Nout', 'likelihoods', 'darFrames', 'partIDs', 'NeiOut', 'pointIndexing', 'pointIDx'
        clear('Vout', 'Nout', 'darFrames', 'NeiOut', 'pointIndexing', 'pointIDx');  % only likelihoods and 'partIDs' are needed
        lenF = size(F,2);        
        fringAll = ComputeFringDeep(F, kkMax, lenF);
        fring = fringAll{kkMax};
        temp = [partIDs(:, 2), likelihoods]';
%         clear('Vout', 'Nout', 'NeiOut', 'pointIndexing', 'pointIDx', 'darFrames', 'NeiOut');
        
        if isempty(crossScaleStructure{kk})
            table = zeros(lenSc, lenF);
        else
            table = crossScaleStructure{kk};
        end
        
        tempF2part = {};
        for ii = 1:size(F,2)
            tempF2part{ii} = [];
        end
        for ii = 1:size(partIDs,1)
            tempF2part{partIDs(ii, 2)} = [tempF2part{partIDs(ii, 2)}, ii];
        end
        
        outLine = zeros(1, lenF);
        facePrev = 0;
        
        for j = 1:lenF

            adjacentFaces = fring{j};
            if isempty(adjacentFaces)
                continue;
            end  
            ids = [tempF2part{adjacentFaces}];
            likeLocal = temp(2, ids);
            if (min(likeLocal) > threshMin) & (sum(likeLocal)/length(likeLocal)>threshAvg)
                outLine(j) = 1;  % no need to proced to larger scales
            end
        end
        
        table(curScale, :) = outLine;
        crossScaleStructure{kk} = table;
        disp(kk);
    end
    
end

