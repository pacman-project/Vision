% This function assigns IDs to each file with statistics
% pairsLeft and pairsRight have the following format:
% part_ID_Central, part_ID_Left, x,y,z, alpha, beta, gamma

function assignIdsToStatisticsFile(pairsLeft, pairsRight, pairClusteringOptions, nPrevClusters, lenFiles, statMapProperties, is_GPU)


    idsD = [(3:8)'; (16:21)'];
    idsPartIds = [1,2,15]';
    idsLeft = [1:6]';
    idsRight = [7:12]';
    
    distThresh = pairClusteringOptions.penaltyThresh * 1.5;
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

    for i = lenFiles:lenFiles %20:20:lenFiles
        
        str = ['Temp/outputStatisticsAll_', num2str(i) ,'.mat'];
        aa = load(str);
        outputStatisticsAll = aa.outputStatisticsAll;
        outputStatistics = [outputStatisticsAll{:}];
        
        % convert to the world's units
        outputStatistics = convertFileToWorldUnits(outputStatistics, xyzStep, angleStep);
        partIds = gpuArray(outputStatistics(idsPartIds, :));
        outputStatistics = outputStatistics(idsD, :);
        
        lenSt = size(outputStatistics, 2);
        outList = zeros(3, lenSt);
        outList(2, :) = gather(partIds(1,:));
        
        for kk = 1:2:3  % left and right pairs
            
            centrals = partIds(1,:);
            if kk == 1 
                pairs = pairsLeft; 
                combs = combsLeft;
                left_right = partIds(2, :);
                idsT = idsLeft;
                clCentrals = clusterCL;
                cl_left_right = clusterLL;
            elseif kk == 3
                pairs = pairsRight; 
                combs = combsRight;
                left_right = partIds(3, :);
                idsT = idsRight;
                clCentrals = clusterCR;
                cl_left_right = clusterRR;
            end
            
            for j = 1:size(combs,1)
                % find this pair in outputStatistics
                ids = centrals == combs(j,1) & left_right == combs(j,2);
                
                % all statistical points 
                outputStatisticsT = outputStatistics(idsT, ids);
                % cluster centres
                idsCl = clCentrals == combs(j,1) & cl_left_right == combs(j,2);
                idsClN = find(idsCl);
                clusterCentres = pairs(idsCl, 3:8);
                D = Mixed_Eucl_Quat_dist(outputStatisticsT, clusterCentres', alpha);
                [M,I] = min(D,[],2);
                if kk == 1
                    outList(kk, ids) = idsClN(I);                
                else
                    outList(kk, ids) = idsClN(I) + numPairsL;
                end
                
                disp(j);
            end
        end
        strOut = ['Temp/outList_',num2str(i), '.mat'];
        save(strOut, 'outList');
    end
end

function combs = uniqueMyRows(matr)
    values = 9999 * matr(1,:) + matr(2,:);
    [~,ia,~] = unique(values);
    combs = matr(:, ia');
end

function statistics = convertFileToWorldUnits(statistics, xyzStep, angleStep)
    idsCoords = [(3:5)'; (16:18)'];
    idsAngles = [(6:8)'; (19:21)'];
    statistics = double(statistics);
    statistics(idsCoords, :) = statistics(idsCoords, :) * xyzStep;
    statistics(idsAngles, :) = statistics(idsAngles, :) * angleStep;
end