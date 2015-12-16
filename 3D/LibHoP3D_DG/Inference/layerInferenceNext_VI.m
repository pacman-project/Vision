% this is to do fast inference from the statistical map

function is_ok = layerInferenceNext_VI(list_input, list_els, lenFiles, layerID, is_GPU_USED)

    is_visualization = 1;

    str = ['Temp/partSelection', num2str(layerID), '.mat'];
    dd = load(str); % 'partsOut', 'coverageOut', 'lenOut', 'pairList';
    numParts = dd.nNClusters{layerID};
    partsOut = dd.partsOut;
    pairsAll = dd.pairsAll;
    clear('dd');

    
    for i = lenFiles:lenFiles  % 20:20:lenFiles
        
        str = ['Temp/outputStatisticsAll_', num2str(i) ,'.mat'];
        aa = load(str);
        outputStatisticsAll = aa.outputStatisticsAll;
        outputStatistics = [outputStatisticsAll{:}];
        
        str = ['Temp/outputCoordsAll_', num2str(i) ,'.mat'];
        aa = load(str);
        outputCoords = aa.outputCoordsAll;
        outputCoords = outputCoords{1};
        
        strIn = ['Temp/outList_',num2str(i), '.mat'];
        aa = load(strIn);
        outList = aa.outList;      
        
        fileName = list_els{i};
        inFileM = [fileName(1:end-4), '.mat'];     % 'Vout', 'Nout', 'darFrames', 'partIDs', 'pointIDx'
        inFilePS = [fileName(1:end-4), 'PS.mat'];  % 'NeiOut'
%         inFileLH = [fileName(1:end-4), 'LH.mat'];  % 'likelihoods', 'partIDs'
%         inFileNP = [fileName(1:end-4), 'NP.mat'];  % 'numPoints'
        
        aaa = load(inFileM);
        Vout = aaa.Vout;
        Nout = aaa.Nout;
        darFrames = aaa.darFrames;
        partIDs = aaa.partIDs;
        pointIDx = aaa.pointIDx;
        
        aaa = load(inFilePS);
        NeiOut = aaa.NeiOut;
        clear('aaa');
        
        if is_GPU_USED
            outList = gpuArray(outList);
            partsCand = zeros(1, size(outList, 2));
        end
        for k = 1:numParts
            ids1 = outList(1, :) == partsOut(k, 1);
            ids2 = outList(2, :) == partsOut(k, 2);
            ids3 = outList(3, :) == partsOut(k, 3);
            ids1 = ids1 & ids2 & ids3;
            partsCand(ids1) = k;
        end
        
        partsCand = gather(partsCand);
        partsNext = zeros(10^6,2);
        partsNeighOutNext = zeros(10^6,3);
        numPartsNext = 0;
        
        centralPartIDx = outputCoords(2, :);
        lenT = length(centralPartIDx);
        k = 1;
        startID = 1;
        prevCentr = centralPartIDx;
        
        while k < lenT
            curCentr = centralPartIDx(k);
            if curCentr ~= prevCentr % do inference here
                endID = k-1;
                pC = partsCand(startID:endID);
                pC = pC(pC > 0);
                numPartsNext = numPartsNext + 1;
                if isempty(pC)
                    partsNext(numPartsNext, 1) = prevCentr;
                    partsNext(numPartsNext, 2) = 0;    % nothing can be inferred
                    partsNeighOutNext(numPartsNext, :) = [0,0,0];
                else
                    p = mode(pC);
                    partsNext(numPartsNext, 1) = prevCentr;
                    partsNext(numPartsNext, 2) = p;
                    idss = find(pC == p);
                    idss = idss(1);
                    idss = idss + startID - 1;
                    partsNeighOutNext(numPartsNext, :) = outputCoords(2:4, idss)';
                end
                
                prevCentr = curCentr;
                startID = k;
            else
                prevCentr = curCentr;
            end
            
            k = k + 1;
            if mod(k, 10^6) == 0
                str = [num2str(k), ' out of ', num2str(lenT)];
                disp(str);
            end
        end
        
        partsNext = partsNext(1:numPartsNext, :);
        ids = partsNext(:, 2) > 0;
        partsNext = partsNext(ids, :);
        partsNeighOutNext = partsNeighOutNext(ids, :);
        outputCoords = outputCoords(:, ids);
        lenPartRealizationsNext = size(partsNext, 1);
        
        tempIdx = zeros(1,size(Vout,2));
        for ii = 1:size(Vout,2)
            tempIdx(ii) = pointIDx{partIDs(ii, 2)}(partIDs(ii, 3));  % pointIDx of all parts in Vout
        end
        
        %% save the output of the files
        
        % visualize the inference results
        strRep =   ['layer', num2str(layerID-1)];
        strRepTo = ['layer', num2str(layerID)];
        
        fileNameNext = strrep(fileName, strRep, strRepTo);
        outFileM =  [fileNameNext(1:end-4), '.mat'];    % 'Vout', 'Nout', 'darFrames', 'partIDs', 'pointIDx'
        outFilePS = [fileNameNext(1:end-4), 'PS.mat'];  % 'NeiOut'
%         outFileLH = [fileNameNext(1:end-4), 'LH.mat'];  % 'likelihoods', 'partIDs'
%         outFileNP = [fileNameNext(1:end-4), 'NP.mat'];  % 'numPoints'

        [~,idxsIntoA] = intersect(tempIdx, partsNext(:,1)');

        Vout = Vout(:,idxsIntoA);
        Nout = Nout(:,idxsIntoA);
        darFrames = darFrames(idxsIntoA, :);
        partIDs = partIDs(idxsIntoA, :);
        partIDs(:, 1) = partsNext(:, 2);
        %  pointIDx remains the same!!!
        
        tic
        NeiOutNext = {};
        if layerID == 5
            for ii = 1:lenPartRealizationsNext
                 bb = unique([NeiOut{partsNeighOutNext(ii, 1)}, NeiOut{partsNeighOutNext(ii, 2)}, NeiOut{partsNeighOutNext(ii, 3)}]);
                 evens = mod(bb,2)==0;
                 NeiOutNext{partsNext(ii, 1)} = bb(evens);
            end
        else
            for ii = 1:lenPartRealizationsNext
                 NeiOutNext{partsNext(ii, 1)} = [NeiOut{partsNeighOutNext(ii, 1)}, NeiOut{partsNeighOutNext(ii, 2)}, NeiOut{partsNeighOutNext(ii, 3)}];
            end
        end

        NeiOut = NeiOutNext;

        if is_visualization

            for kk = 1:numParts
                figure;
                iids = partIDs(:, 1) == kk;
                VV = Vout(:, iids);
                scatter3(VV(1,:), VV(2,:), VV(3,:)');
                a = 2;
            end

        end

        save(outFileM,  'Vout', 'Nout', 'darFrames', 'partIDs', 'pointIDx');    % needed for statistics collection
        save(outFilePS, 'NeiOut');                                              % needed for part selection
%             save(outFileLH, 'likelihoods', 'partIDs');                              % needed for  reviseFaces
%             save(outFileNP, 'numPoints');         

        a = 2;

        
    end


    is_ok = true;
end







