% This function assigns IDs to each file with statistics
% pairsLeft and pairsRight have the following format:
% part_ID_Central, part_ID_Left, x,y,z, alpha, beta, gamma

function assignIdsToStatisticsFile(pairsLeft, pairsRight, pairClusteringOptions, nPrevClusters, lenFiles, statMapProperties, layerID, is_GPU)


    idsD = [(3:9)']; %; (16:21)'];
    idsPartIds = [1;2];
    idsLeftRight = [1:7]';
    
    if layerID <= 5
        numCombs = 6;
    else
        numCombs = 10;
    end
    
    distThresh = 6;
    alpha = pairClusteringOptions.alpha;
    
    numPairsL = size(pairsLeft, 1);

    xyzStep = statMapProperties.xyzStep;
    angleStep = statMapProperties.angleStep; 
    
    % find all pairs of parts that belong to some clusters
    clusterCL = pairsLeft(:,1);
    clusterLL = pairsLeft(:,2);
    combsLeft = unique([clusterCL, clusterLL], 'rows');
    
    clusterCR = pairsRight(:,1);
    clusterRR = pairsRight(:,2);
    combsRight = unique([clusterCR, clusterRR], 'rows');

    parfor i = 1:lenFiles
        
        str = ['Temp/outputStatisticsAll_', num2str(i) ,'.mat'];
        
        if ~exist(str, 'file')
            continue;
        end
        aa = load(str); 
        outputStatistics = aa.outputStatisticsAll;
        
        str = ['Temp/outputCoordsAll_', num2str(i) ,'.mat'];
        aa = load(str);
        outputCoordsAll = aa.outputCoordsAll;
        leftIds  = outputCoordsAll(3, :) > 0;
        rightIds = outputCoordsAll(4, :) > 0;

        
        % it becomes of the following format:partID, partID, x, y, z, Q1, Q2, Q3, Q4
        outputStatistics = double(outputStatistics);
        outputStatistics(6:9, :) = convertToQuaternions(outputStatistics(6:end, :), statMapProperties, layerID);  % 9 values are converted to 1 quaternion value
        outputStatistics = outputStatistics(1:9, :);
        
        % convert to the world's units
        outputStatistics = convertFileToWorldUnits(outputStatistics, xyzStep, angleStep);
        partIds = outputStatistics(idsPartIds, :);
        outputStatistics = outputStatistics(idsD, :);
        
        lenSt = size(outputStatistics, 2);
        outList = zeros(3, lenSt);
        outList(2, :) = gather(partIds(1,:));
        
        left_right = partIds(2, :);
        
        for kk = 1:2:3  % left and right pairs
            
            centrals = partIds(1,:);
            if kk == 1 
                pairs = pairsLeft; 
                combs = combsLeft;
                clCentrals = clusterCL;
                cl_left_right = clusterLL;
                idsCur = leftIds;
            elseif kk == 3
                pairs = pairsRight; 
                combs = combsRight;
                clCentrals = clusterCR;
                cl_left_right = clusterRR;
                idsCur = rightIds;
            end
            
            for j = 1:size(combs,1)
                % find this pair in outputStatistics
                ids = centrals == combs(j,1) & left_right == combs(j,2);
                idsAdd = ids & idsCur;
                % all statistical points 
                outputStatisticsT = outputStatistics(idsLeftRight, idsAdd);
                % cluster centres
                idsCl = clCentrals == combs(j,1) & cl_left_right == combs(j,2);
                idsClN = find(idsCl);
                clusterDescription = pairs(idsCl, 3:34);
%                 D = Mixed_Eucl_Quat_dist(outputStatisticsT, clusterCentres', alpha);
                D = Mixed_Eucl_Quat_cluster_dist(outputStatisticsT', clusterDescription, alpha);
                [M,I] = min(D,[],2);
                idsAdd1 = find(idsAdd);
                idsAdd2 = find(M <= distThresh);
%                 outputCoordsAllT(idsLeftRight, idsAdd) = ;
%                 
%                 if kk == 1
%                     outList(kk, idsAdd) = idsClN(I);                
%                 else
%                     outList(kk, idsAdd) = idsClN(I) + numPairsL;
%                 end
                

                if kk == 1
                    outList(kk, idsAdd1(idsAdd2)) = idsClN(I(idsAdd2));
%                     outputCoordsAllT(2:3, idsAdd1(idsAdd2)) = outputCoordsAll(2:3, idsAdd1(idsAdd2));
                else
                    outList(kk, idsAdd1(idsAdd2)) = idsClN(I(idsAdd2)) + numPairsL;
%                     idsRR = [2,4]';
%                     outputCoordsAllT(idsRR, idsAdd1(idsAdd2)) = outputCoordsAll(idsRR, idsAdd1(idsAdd2));
                end
                
%                 numPP = length(idsAdd1(idsAdd2));
%                 
%                 % add the second hypothesis (if it is less than the threshold)
%                 rows = (1:length(I))';
%                 indsInD = sub2ind(size(D), rows, I);
%                 D(indsInD) = 100;
%                 [M,I] = min(D,[],2);
%                 idsAdd2 = find(M <= distThresh);
%                 if kk == 1
%                     outList(kk, numPP + idsAdd1(idsAdd2)) = idsClN(I(idsAdd2));
%                     outputCoordsAllT(2:3, numPP + idsAdd1(idsAdd2)) = outputCoordsAll(2:3, idsAdd1(idsAdd2));
%                 else
%                     outList(kk, numPP + idsAdd1(idsAdd2)) = idsClN(I(idsAdd2)) + numPairsL;
%                     outputCoordsAllT(idsRR, numPP + idsAdd1(idsAdd2)) = outputCoordsAll(idsRR, idsAdd1(idsAdd2));
%                 end
            end
        end
        outputCoordsAllT = outputCoordsAll;
        
        % re-arrange outList to make all triples
%         tt = [outList; outputCoordsAllT(2:4, :)];
        
        lenSt = size(outList, 2);
        idxPrev = outputCoordsAllT(2,1);
        PIcentralPrev = outList(2,1);
        beginPrev = 1;
        outListTemp = {};
        outCoordsTemp = {};
        
        disp('Making triples');
        for ii  = 2:lenSt
            if idxPrev ~= outputCoordsAllT(2,ii)
                endPrev = ii - 1;
                leftPIs =  outList(1, beginPrev:endPrev); %outputCoordsAll(3, beginPrev:endPrev);
                rightPIs = outList(3, beginPrev:endPrev); %outputCoordsAll(4,beginPrev:endPrev);
                leftIDs = find(leftPIs);
                rightIDs = find(rightPIs);
                leftPIs = leftPIs(leftIDs);
                rightPIs = rightPIs(rightIDs);
                leftOC = outputCoordsAllT(3, leftIDs + beginPrev-1);
                rightOC =outputCoordsAllT(4, rightIDs + beginPrev-1);
                
                numLeft = length(leftPIs);
                numRight = length(rightPIs);
                if numLeft > numCombs
                    idsCurRF = randperm(numLeft, numCombs);
                    leftPIs= leftPIs(idsCurRF);
                    leftOC = leftOC(idsCurRF);
                    numLeft = numCombs;
                end
                if numRight > numCombs
                    idsCurRF = randperm(numRight, numCombs);
                    rightPIs= rightPIs(idsCurRF);
                    rightOC = rightOC(idsCurRF);
                    numRight = numCombs;
                end
                
                vec1 = 1:numLeft; vec2 = 1:numRight;
                [p,q] = meshgrid(vec1, vec2);
                vec1 = p(:);
                vec2 = q(:);
                lenVec = length(vec1);
                olTemp = int32(ones(3, lenVec));  % out list
                ocTemp = int32(ones(4, lenVec));  % output coords
                
                olTemp(2, :) = olTemp(2, :) * PIcentralPrev;
                olTemp(1, :) = leftPIs(vec1)';
                olTemp(3, :) = rightPIs(vec2)';
                ocTemp(2, :) = ocTemp(2, :) * idxPrev;
                ocTemp(3, :) = leftOC(vec1)';
                ocTemp(4, :) = rightOC(vec2)';
                
                outListTemp{ii} = olTemp;
                outCoordsTemp{ii} = ocTemp;

                beginPrev = ii;
                idxPrev = outputCoordsAllT(2,ii);
                PIcentralPrev = outList(2,ii);
            end
        end

        outList = int16([outListTemp{:}]);
        outputCoordsAll_A = [outCoordsTemp{:}];
        a = 2;
        
        str = ['File ', num2str(i)];
        disp(str);
        strOut1 = ['Temp/outList_',num2str(i), '.mat'];
        parSaveOL(strOut1, outList, 'outList');
        
        strOut2 = ['Temp/outputCoordsAll_A', num2str(i) ,'.mat'];
        parSaveOC(strOut2, outputCoordsAll_A, 'outputCoordsAll_A');
    end
end

function parSaveOL(fileName, outList, varName)
    save(fileName, varName);
end

function parSaveOC(fileName, outputCoordsAll_A, varName)
    save(fileName, varName);
end

function combs = uniqueMyRows(matr)
    values = 9999 * matr(1,:) + matr(2,:);
    [~,ia,~] = unique(values);
    combs = matr(:, ia');
end

% function combs = uniqueMyRows3(matr)
%     values = 99999999 * matr(1,:) + 9999 * matr(2,:) + matr(3,:);
%     [~,ia,~] = unique(values);
%     combs = matr(:, ia');
% end

% convert normals or axis to quaternions
% convert normals or frames to quaternions
function statisticsOut = convertToQuaternions(statistics, statMapProperties, layerID)
    angleStep = statMapProperties.angleStep;
    statistics = statistics * angleStep;
    
    if layerID < 5 % convert NORMAL to quaternions 
        normsV = sqrt(statistics(1,:).^2 + statistics(2,:).^2 + statistics(3,:).^2); % normalize the vectors
        statistics(1:3,:) = statistics(1:3,:)./repmat(normsV, [3,1]);
        % compute quaternion for each vector
        statisticsOut = computeVectorQuaternionVectorized(statistics(1:3,:), [0,0,1]);
        statisticsOut = qconj(statisticsOut);
        
    else  % convert FRAME to quaternion
        normsV = sqrt(statistics(1,:).^2 + statistics(2,:).^2 + statistics(3,:).^2);
        statistics(1:3,:) = statistics(1:3,:)./repmat(normsV, [3,1]);
        normsX = sqrt(statistics(4,:).^2 + statistics(5,:).^2 + statistics(6,:).^2);
        statistics(4:6,:) = statistics(4:6,:)./repmat(normsX, [3,1]);
        normsY = sqrt(statistics(7,:).^2 + statistics(8,:).^2 + statistics(9,:).^2);
        statistics(7:9,:) = statistics(7:9,:)./repmat(normsY, [3,1]);
        
        tic
        lenStat = size(statistics, 2);
        statisticsOut = zeros(4, lenStat);
            parfor i = 1:lenStat
                dcm = [statistics(4:6,i), statistics(7:9,i), statistics(1:3,i)]'; % [X; Y; N]  
                q = dcm2q(dcm);
                statisticsOut(:, i) = q;
            end
            statisticsOut = qnorm(statisticsOut);
        toc
    end
end
% function statisticsOut = convertToQuaternions(statistics, statMapProperties, layerID)
%     angleStep = statMapProperties.angleStep;
%     statistics = statistics * angleStep;
%     
%     if layerID < 6 % convert Normal to quaternions 
%         normsV = sqrt(statistics(1,:).^2 + statistics(2,:).^2 + statistics(3,:).^2); % normalize the vectors
%         statistics(1:3,:) = statistics(1:3,:)./repmat(normsV, [3,1]);
%         
%         % compute quaternion for each vector
%         statisticsOut = computeVectorQuaternionVectorized(statistics(1:3,:), [0,0,1]);
%         statisticsOut = qconj(statisticsOut);
%     else  % convert frame to quaternion
%         normsV = sqrt(statistics(1,:).^2 + statistics(2,:).^2 + statistics(3,:).^2);
%         statistics(1:3,:) = statistics(1:3,:)./repmat(normsV, [3,1]);
%         normsV = sqrt(statistics(4,:).^2 + statistics(5,:).^2 + statistics(6,:).^2);
%         statistics(4:6,:) = statistics(1:3,:)./repmat(normsV, [3,1]);
%         normsV = sqrt(statistics(7,:).^2 + statistics(8,:).^2 + statistics(9,:).^2);
%         statistics(7:9,:) = statistics(7:9,:)./repmat(normsV, [3,1]);
%     end
%     a = 2;
% end


function statistics = convertFileToWorldUnits(statistics, xyzStep, angleStep)
    idsCoords = [(3:5)']; % (16:18)'];
    idsAngles = [(6:9)']; % (19:21)'];
    statistics = double(statistics);
    statistics(idsCoords, :) = statistics(idsCoords, :) * xyzStep;
%     statistics(idsAngles, :) = statistics(idsAngles, :) * angleStep;
end