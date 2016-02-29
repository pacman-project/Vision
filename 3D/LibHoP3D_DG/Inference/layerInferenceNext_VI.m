% this is to do fast inference from the statistical map

function is_ok = layerInferenceNext_VI(list_input, list_els, lenFiles, layerID, is_GPU_USED)

    is_visualization = 1;

    % str = ['Temp/partSelection', num2str(layerID), '.mat'];
    str = ['Temp/','layer', num2str(layerID), '/partSelection', num2str(layerID), '.mat']
    dd = load(str); % 'partsOut', 'coverageOut', 'lenOut', 'pairList';
    numParts = dd.nNClusters{layerID};
    partsOut = dd.partsOut;
    pairsAll = dd.pairsAll;
    clear('dd');
    

    strLayer = ['Temp/OR_node_layer_', num2str(layerID)];
    dd = load(strLayer);
    ORTable = dd.ORTable; 

    
    numPartsBeforeORNodes = size(partsOut, 1); 
    if layerID == 4 || layerID == 6
        numParts = numParts - 1; % first do inference of those parts that are not 
    end
    
    for i = 1:1%lenFiles
        
%         str = ['Temp/outputStatisticsAll_', num2str(i) ,'.mat'];
%         aa = load(str);
%         outputStatisticsAll = aa.outputStatisticsAll;
        
        str = ['Temp/outputCoordsAll_A', num2str(i) ,'.mat'];
        if ~exist(str, 'file')
            continue;
        end
        aa = load(str);
        outputCoords = aa.outputCoordsAll_A;
        
        strIn = ['Temp/outList_',num2str(i), '.mat'];
        aa = load(strIn);
        outList = aa.outList;      
        
        if layerID == 3
            sc = 1;
            strScale = ['scale_', num2str(sc)];
        else
            strScale = [];
        end
        
        fileName = list_els{i};
        inFileM = [fileName(1:end-4),strScale, '.mat'];     % 'Vout', 'Nout', 'darFrames', 'partIDs', 'pointIDx'
        inFilePS = [fileName(1:end-4), strScale, 'PS.mat'];  % 'NeiOut'
%         inFileLH = [fileName(1:end-4), 'LH.mat'];  % 'likelihoods', 'partIDs'
%         inFileNP = [fileName(1:end-4), 'NP.mat'];  % 'numPoints'
        
        aaa = load(inFileM);
        Vout = aaa.Vout;
        Nout = aaa.Nout;
        darFrames = aaa.darFrames;
        partIDs = aaa.partIDs;
        partIDs(:, 5) = (1:size(partIDs, 1))';
        
%         aaa = load(inFilePS);
%         NeiOut = aaa.NeiOut;
%         clear('aaa');
        
        if is_GPU_USED
            outList = gpuArray(outList);
            partsCand = zeros(1, size(outList, 2));
        end
        for k = 1:numPartsBeforeORNodes
            ids1 = outList(1, :) == partsOut(k, 1);
            ids2 = outList(2, :) == partsOut(k, 2);
            ids3 = outList(3, :) == partsOut(k, 3);
            ids1 = ids1 & ids2 & ids3;
            partsCand(ids1) = k;
        end
        
        partsCand = gather(partsCand);
        partsNext = zeros(10^6,2);
%         partsNeighOutNext = zeros(10^6,3);
        numPartsNext = 0;
        
        centralRealizationID = partIDs(outputCoords(2, :),5);
        lenT = length(centralRealizationID);
        k = 1;
        startID = 1;
        prevCentrRelaizID = centralRealizationID;
        
%         while k < lenT
%             curCentr = centralRealizationID(k);
%             if curCentr ~= prevCentrRelaizID % do inference here
%                 endID = k-1;
%                 pC = partsCand(startID:endID);
%                 pC = pC(pC > 0);
%                 numPartsNext = numPartsNext + 1;
%                 if isempty(pC)
%                     partsNext(numPartsNext, 1) = prevCentrRelaizID;
%                     partsNext(numPartsNext, 2) = 0;    % nothing can be inferred
%                     partsNeighOutNext(numPartsNext, :) = [0,0,0];
%                 else
%                     p = mode(pC);
%                     partsNext(numPartsNext, 1) = prevCentrRelaizID;
%                     partsNext(numPartsNext, 2) = p;
%                     idss = find(pC == p);
%                     idss = idss(1);
%                     idss = idss + startID - 1;
%                     partsNeighOutNext(numPartsNext, :) = outputCoords(2:4, idss)';
%                 end
%                 
%                 prevCentrRelaizID = curCentr;
%                 startID = k;
%             else
%                 prevCentrRelaizID = curCentr;
%             end
%             
%             k = k + 1;
%             if mod(k, 10^6) == 0
%                 str = [num2str(k), ' out of ', num2str(lenT)];
%                 disp(str);
%             end
%         end
        
        partsNext = partsNext(1:numPartsNext, :);
        ids = partsNext(:, 2) > 0;
        partsNext = partsNext(ids, :);
%         partsNeighOutNext = partsNeighOutNext(ids, :);
%         outputCoords = outputCoords(:, ids);
        lenPartRealizationsNext = size(partsNext, 1);
           

        partsNext(:, 2) = ORTable(partsNext(:, 2));
        
        %% save the output of the files
        
        % visualize the inference results
        strRep =   ['layer', num2str(layerID-1)];
        strRepTo = ['layer', num2str(layerID)];
        
        fileNameNext = strrep(fileName, strRep, strRepTo);
        outFileM =  [fileNameNext(1:end-4), '.mat'];    % 'Vout', 'Nout', 'darFrames', 'partIDs', 'pointIDx'
        outFilePS = [fileNameNext(1:end-4), 'PS.mat'];  % 'NeiOut'
%         outFileLH = [fileNameNext(1:end-4), 'LH.mat'];  % 'likelihoods', 'partIDs'
%         outFileNP = [fileNameNext(1:end-4), 'NP.mat'];  % 'numPoints'

        tempRealizationIDs = partIDs(:, 5);
        [~,idxsIntoA] = intersect(tempRealizationIDs, partsNext(:,1)');
        
        Vout = Vout(:,idxsIntoA);
        Nout = Nout(:,idxsIntoA);
        darFrames = darFrames(idxsIntoA, :);
        partIDs = partIDs(idxsIntoA, :);
        
        if layerID == 4 || layerID == 6
            partsNext(:, 2) = partsNext(:, 2) + 1;
        end
        partIDs(:, 1) = partsNext(:, 2);
        %  pointIDx remains the same!!!
        

%         for ii = 1:lenPartRealizationsNext
%              NeiOutNext{ii} = [NeiOut{partsNeighOutNext(ii, 1)}, NeiOut{partsNeighOutNext(ii, 2)}, NeiOut{partsNeighOutNext(ii, 3)}];
%         end

        
        if layerID == 4  % deleting points with scale == 1
            sc = 2;
            strScale =  ['scale_', num2str(sc)];
            fileName =  list_els{i};
            inFileM1 =  [fileName(1:end-4), strScale, '.mat'];   % 'Vout', 'Nout', 'darFrames', 'partIDs', 'pointIDx'
            inFilePS1 = [fileName(1:end-4), strScale, 'PS.mat']; % 'NeiOut'
            inFileNP1 =  [fileName(1:end-4), 'scale_1NP.mat'];
            
            inFileM1 = strrep(inFileM1, 'layer3', 'layer2');
            inFilePS1 = strrep(inFilePS1, 'layer3', 'layer2');
            inFileNP1 = strrep(inFileNP1, 'layer3', 'layer2');
            
            
            tt = load(inFileNP1);
            aa = load(inFileM1);
            Vout = [Vout, aa.Vout];
            Nout = [Nout, aa.Nout];
            darFrames = [darFrames; aa.darFrames];
            partIDs = [partIDs(:, 1:4); aa.partIDs];
            
            aa = load(inFilePS1);
            NO = aa.NeiOut;
            pointScales = tt.pointScales;
            
            
%             for ii = 1:size(Vout, 2)
%                 if ii <= lenPartRealizationsNext
%                     temp = NeiOutNext{ii};
%                     scales = pointScales(temp);
%                     scales = scales > 1;
%                     temp = temp(scales);
%                     NeiOut{ii} = temp;
%                 else
%                     NeiOut{ii} = NO{ii - lenPartRealizationsNext};
%                 end
%             end  
        elseif layerID == 6
            sc = 3;
            strScale =  ['scale_', num2str(sc)];
            fileName =  list_els{i};
            inFileM1 =  [fileName(1:end-4), strScale, '.mat'];   % 'Vout', 'Nout', 'darFrames', 'partIDs', 'pointIDx'
            inFilePS1 = [fileName(1:end-4), strScale, 'PS.mat']; % 'NeiOut'
            inFileNP1 =  [fileName(1:end-4), 'scale_1NP.mat'];
            
            inFileM1 = strrep(inFileM1,   'layer5', 'layer2');
            inFilePS1 = strrep(inFilePS1, 'layer5', 'layer2');
            inFileNP1 = strrep(inFileNP1, 'layer5', 'layer2');
            
            tt = load(inFileNP1);
            aa = load(inFileM1);
            Vout = [Vout, aa.Vout];
            Nout = [Nout, aa.Nout];
            darFrames = [darFrames; aa.darFrames];
            partIDs = [partIDs(:, 1:4); aa.partIDs];
            
            aa = load(inFilePS1);
            NO = aa.NeiOut;
            pointScales = tt.pointScales;
            
%             for ii = 1:size(Vout, 2)
%                 if ii <= lenPartRealizationsNext
%                     temp = NeiOutNext{ii};
%                     scales = pointScales(temp);
%                     scales = scales > 2;
%                     temp = temp(scales);
%                     NeiOut{ii} = temp;
%                 else
%                     NeiOut{ii} = NO{ii - lenPartRealizationsNext};
%                 end
%             end  
        else
%             NeiOut = NeiOutNext;
        end

        if is_visualization
            
            vecLen = 0.1;
            vectColors = eye(3);
            
            if layerID == 4 || layerID == 6
                numPartsOut = numParts + 1; % first do inference of those parts that are not 
            else
                numPartsOut = numParts;
            end
            for kk = 1:numPartsOut
                figure;
                iids = find(partIDs(:, 1) == kk);
                VV = Vout(:, iids);
                scatter3(VV(1,:), VV(2,:), VV(3,:)');

%                 hold on
%                 lenParts = length(iids);
%                 numNormals = min(100, lenParts);
%                 rr = randperm(lenParts, numNormals);
% 
%                 for kkkk = 1:numNormals
%                     VD = darFrames(rr(kkkk), :);
%                     VD = reshape(VD, [3,3]);
%                     Vcentre = Vout(:, rr(kkkk));
%                     plotFrame(VD, vecLen, vectColors, Vcentre);
%                     a = 2;
%                 end
%                 
%                 a = 2;
            end

        end
        
        % check if the folder exists
        [folderName, ~] = getFolderName(outFileM);
        
        if ~exist(folderName, 'dir')
           mkdir(folderName);
        end

%         save(outFileM,  'Vout', 'Nout', 'darFrames', 'partIDs');    % needed for statistics collection
        parSave(outFileM, 'Vout', 'Nout', 'darFrames', 'partIDs', Vout, Nout, darFrames, partIDs);

%         save(outFilePS, 'NeiOut');                                              % needed for part selection
%             save(outFileLH, 'likelihoods', 'partIDs');                              % needed for  reviseFaces
%             save(outFileNP, 'numPoints');         

        a = 2; 
    end


    is_ok = true;
end

function parSave(outFileM, str1, str2, str3, str4, Vout, Nout, darFrames, partIDs)
    save(outFileM, str1, str2, str3, str4);
end

function plotFrame(V, vecLen, vectColors, point)
    for ii = 1:3
        curVect = V(:, ii);
        curColor = vectColors(ii, :);
        XX = [point(1), point(1) + vecLen * curVect(1)];
        YY = [point(2), point(2) + vecLen * curVect(2)];
        ZZ = [point(3), point(3) + vecLen * curVect(3)];

        plot3(XX, YY, ZZ, 'Color', curColor);
        hold on
    end
end


            
%             counter1 = 0;
%             counter2 = 0;
%             counter3 = 0;
%             for ii = 1:lenPartRealizationsNext
%                 if ~isempty(NeiOutNext{ii}) && ~isempty(NO{ii})
%                     NeiOut{ii} = NO{ii};
%                     counter1 = counter1 + 1;
%                 else
%                     if ~isempty(NO{ii})
%                         NeiOut{ii} = NO{ii};
%                         counter2 = counter2 + 1;
%                     end
%                     if ~isempty(NeiOutNext{ii})
%                         % only the points with pointScales > 1 should survive
%                         temp = NeiOutNext{ii};
%                         scales = pointScales(temp);
%                         scales = scales > 1;
%                         temp = temp(scales);
%                         NeiOut{ii} = temp;
%                         counter3 = counter3 + 1;
%                     end
%                 end
%             end







