% this function builds statistical maps out of co-occurrence statistic
% files

function [statMapLeft, statMapRight] = buildStatMap_NextLayers(fileIds, fileStart, lenFiles, nPrevClusters, offsetConventional, statMapProperties, layerID, sieveThresh, is_GPU)
                                                    
    idsLeftRight  = (3:14)';
    idsLeftRightQ  = (3:9)';
%     idsXYLeft = (3:4)';
%     idsZLeft = 5;
%     idsRight = (16:21)';
%     idsXYRight = (16:17)';
%     idsZRight = 18;
    
    statMapLeft = {};
    statMapRight = {};
    
    for i = 1:nPrevClusters
        for j = 1:nPrevClusters
            statMapLeft{i,j} = {};
            statMapRight{i,j} = {};
        end
    end
    
    if fileStart ~= 1
        aa = load('Temp/statMap.mat');
        statMapLeft = aa.statMapLeft;
        statMapRight = aa.statMapRight;
    end
    
    [sizeXY, sizeZ, sizeAngle, centreXY, centreZ, centreAngle] = computeStatMapSizes(statMapProperties, offsetConventional);
%     mapSize = sizeXY*sizeXY*sizeZ*sizeAngle*sizeAngle*sizeAngle*sizeAngle;
    
    
    for i = fileStart:lenFiles
        
        str = ['Temp/outputStatisticsAll_', num2str(fileIds(i)) ,'.mat'];
        
        if ~exist(str)
            continue;
        end
        aa = load(str);
        outputStatistics = aa.outputStatisticsAll;
        
        str = ['Temp/outputCoordsAll_', num2str(fileIds(i)) ,'.mat'];
        aa = load(str);
        outputCoordsAll = aa.outputCoordsAll;

        
        if isempty(outputStatistics)
            continue;
        end

        if ~is_GPU
            centrals = outputStatistics(1,:);
            leftIds  = outputCoordsAll(3, :) > 0;
            rightIds = outputCoordsAll(4, :) > 0;
            lefts = outputStatistics(2, leftIds);
            rights = outputStatistics(2, rightIds);
            combsLeft = unique([centrals(leftIds); lefts]', 'rows');
            combsRight = unique([centrals(rightIds); rights]', 'rows');  %TOO SLOW
            
            disp('ERROR: Computations are done on CPU. No filtering by SieveThresh1 is applied');
        else 
            leftIds  = outputCoordsAll(3, :) > 0;
            rightIds = outputCoordsAll(4, :) > 0;
            centrals = gpuArray(outputStatistics(1,:));
            lefts = gpuArray(outputStatistics(2, leftIds));
            rights = gpuArray(outputStatistics(2, rightIds));
            combsLeft = gather(uniqueMyRows(double( [centrals(leftIds); lefts]),   sieveThresh))';
            combsRight = gather(uniqueMyRows(double([centrals(rightIds); rights]), sieveThresh))';
        end
 
        % convert outputStatistics to the units of the statistical map
        
        % it becomes of the following format:partID, partID, x, y, z, Q1, Q2, Q3, Q4
        outputStatistics = double(outputStatistics);
        outputStatistics(6:9, :) = convertToQuaternions(outputStatistics(6:end, :), statMapProperties, layerID);  % 9 values are converted to 1 quaternion value
        outputStatistics = outputStatistics(1:9, :);
        outputStatistics = fromFileToStatMap(outputStatistics, statMapProperties, centreXY, centreZ, centreAngle, offsetConventional, leftIds, rightIds, layerID);

        for k = 1:size(combsLeft, 1)
            if isempty(statMapLeft{combsLeft(k,1), combsLeft(k,2)})
                statMapLeft{combsLeft(k,1), combsLeft(k,2), 1} = [];  % ids
                statMapLeft{combsLeft(k,1), combsLeft(k,2), 2} = [];  % frequencies
            end
        end
        for k = 1:size(combsRight, 1)
            if isempty(statMapRight{combsRight(k,1), combsRight(k,2)})
                statMapRight{combsRight(k,1), combsRight(k,2), 1} = []; % ids
                statMapRight{combsRight(k,1), combsRight(k,2), 2} = []; % frequencies
            end
        end

        
        leftRights = outputStatistics(2,:);
        for jj = 1:2 % left and right combs
            
            if jj == 1
                combsCur = combsLeft;
                statMapCur = statMapLeft;
                idsCur = leftIds;
            else
                combsCur = combsRight;
                statMapCur = statMapRight;
                idsCur = rightIds;
            end
            
            for k = 1:size(combsCur,1)
                disp(k);
                ids = centrals == combsCur(k,1) & leftRights == combsCur(k,2) & idsCur;
                
                statisticsTemp = double(outputStatistics(idsLeftRightQ, ids));
                statisticsTemp = CheckStatistics(statisticsTemp, sizeXY, sizeZ, sizeAngle);

                inds = sub2ind([sizeXY,sizeXY,sizeZ,sizeAngle,sizeAngle,sizeAngle,sizeAngle], statisticsTemp(1,:), statisticsTemp(2,:), statisticsTemp(3,:) ,...
                    statisticsTemp(4,:), statisticsTemp(5,:), statisticsTemp(6,:), statisticsTemp(7,:));
                t = unique(inds);
                if length(t) > 1
                    [freq, idds] = hist(inds, t);
                else
                    idds = t;
                    freq = length(inds);
                end
                if ~isempty(statMapCur{combsCur(k,1), combsCur(k,2), 1}) 
                    iddsOld =  statMapCur{combsCur(k,1), combsCur(k,2), 1};
                    freqOld =  statMapCur{combsCur(k,1), combsCur(k,2), 2};
%                     tempMatr(iddsOld) = freqOld;
%                     tempMatr(idds) = tempMatr(idds) + freq;
%                     idds = find(tempMatr);
%                     freq = tempMatr(idds);

                    idssTemp =  [iddsOld, idds];
                    freqTemp =  [freqOld, freq];
                    [idds, ~, ic] = unique(idssTemp);
                    freq = accumarray(ic, freqTemp);
                    if size(freq, 1) > 1
                        freq = freq';
                    end
                end
                statMapCur{combsCur(k,1), combsCur(k,2), 1} = idds;
                statMapCur{combsCur(k,1), combsCur(k,2), 2} = freq;
            end
            
            if jj == 1
                statMapLeft = statMapCur;
            else
                statMapRight = statMapCur;
            end
        end
        str = ['File ', num2str(i), 'is finished'];
        disp(str);
    end
    
end

function combs = uniqueMyRows(matr, thresh)
    values = 9999 * matr(1,:) + matr(2,:);
    [valuesUnique, ia, ~] = unique(values);
%     [iaSort, ids] = sort(ia);
    if length(valuesUnique) > 1
        freq = hist(values, valuesUnique);
    elseif length(valuesUnique) == 1
        freq = length(values);
    end
    ids = freq>thresh;
    ia = ia(ids);
    combs = matr(:, ia');
end

% unit conversion
function outputStatistics = fromFileToStatMap(outputStatistics, statMapProperties, centreXY, centreZ, centreAngle, offsetConventional, idsLeft, idsRight, layerID)

    idsAngles = (6:9)'; % quaternions
    
    xyzStep = statMapProperties.xyzStep;
    quaternionSteps = statMapProperties.quaternionSteps;
    
    if layerID == 3 || layerID == 5 || layerID == 7 || layerID == 9
        outputStatistics(idsAngles, :) = round(outputStatistics(idsAngles, :)/quaternionSteps) + centreAngle;
        outputStatistics(3, idsLeft)  = outputStatistics(3, idsLeft) - offsetConventional/xyzStep + centreXY;
        outputStatistics(4, idsLeft)  = outputStatistics(4, idsLeft) + centreXY;
        outputStatistics(5, idsLeft)  = outputStatistics(5, idsLeft) + centreZ;
        outputStatistics(3, idsRight) = outputStatistics(3, idsRight) + offsetConventional/xyzStep + centreXY;
        outputStatistics(4, idsRight) = outputStatistics(4, idsRight) + centreXY;
        outputStatistics(5, idsRight) = outputStatistics(5, idsRight) + centreZ;
        
    elseif layerID == 4 || layerID == 6 || layerID == 8 || layerID == 10
        outputStatistics(idsAngles, :) = round(outputStatistics(idsAngles, :)/quaternionSteps) + centreAngle;
        outputStatistics(3, idsLeft)  = outputStatistics(3, idsLeft) + centreXY;
        outputStatistics(4, idsLeft)  = outputStatistics(4, idsLeft) - offsetConventional/xyzStep + centreXY;
        outputStatistics(5, idsLeft)  = outputStatistics(5, idsLeft) + centreZ;
        outputStatistics(3, idsRight) = outputStatistics(3, idsRight) + centreXY;
        outputStatistics(4, idsRight) = outputStatistics(4, idsRight) + offsetConventional/xyzStep + centreXY;
        outputStatistics(5, idsRight) = outputStatistics(5, idsRight) + centreZ;
        
        a = 2;
    end
    
    outputStatistics = round(outputStatistics);
end

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
            % normalize all quaternions
            statisticsOut = qnorm(statisticsOut);   
        toc
        
        a = 2;
    end
end



function statisticsTemp = CheckStatistics(statisticsTemp, sizeXY, sizeZ, sizeAngle)

    idsWrong =             statisticsTemp(1,:) > sizeXY    | statisticsTemp(1,:) < 1;
    idsWrong = idsWrong | (statisticsTemp(2,:) > sizeXY    | statisticsTemp(2,:) < 1);
    idsWrong = idsWrong | (statisticsTemp(3,:) > sizeZ     | statisticsTemp(3,:) < 1);
    idsWrong = idsWrong | (statisticsTemp(4,:) > sizeAngle | statisticsTemp(4,:) < 1);
    idsWrong = idsWrong | (statisticsTemp(5,:) > sizeAngle | statisticsTemp(5,:) < 1);
    idsWrong = idsWrong | (statisticsTemp(6,:) > sizeAngle | statisticsTemp(6,:) < 1);
    idsWrong = idsWrong | (statisticsTemp(7,:) > sizeAngle | statisticsTemp(7,:) < 1);
    
    statisticsTemp = statisticsTemp(:, ~idsWrong);
end




%     % filtering of the statistical maps
%     for i = 1:nPrevClusters
%         for j = 1:nPrevClusters 
%             if ~isempty(statMapLeft{i,j,1})
%                 ids = statMapLeft{i,j,1};
%                 freqs = statMapLeft{i,j,2};
%                 
%                 sizeMatr = [sizeXY,sizeXY,sizeZ,sizeAngle,sizeAngle,sizeAngle];
%                 
%                 if length(ids) > 1
%                     [i1,i2,i3,i4,i5,i6] = ind2sub(sizeMatr, idds);
%                     points = [i1;i2;i3;i4;i5;i6]';
%                     
%                 end
%             end
%             
%             if ~isempty(statMapRight{i,j,1})
%                 ids = statMapRight{i,j,1};
%                 freqs = statMapRight{i,j,1};
%                 
% 
%                 if length(ids) == 1
%                     continue;
%                 end
%             end
%         end
%     end
%     
%     
%     statMapLeft, statMapRight












