% this function builds statistical maps out of co-occurrence statistic
% files

function [statMapLeft, statMapRight] = buildStatMap_NextLayers(lenFiles, nPrevClusters, offsetConventional, statMapProperties, layerID, is_GPU)
                                                    
    idsLeft  = (3:8)';
%     idsXYLeft = (3:4)';
%     idsZLeft = 5;
    idsRight = (16:21)';
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
    
    [sizeXY, sizeZ, sizeAngle, centreXY, centreZ, centreAngle] = computeStatMapSizes(statMapProperties, offsetConventional);
    
    
    for i = lenFiles:lenFiles %20:20:lenFiles
        
        str = ['Temp/outputStatisticsAll_', num2str(i) ,'.mat'];
        aa = load(str);
        outputStatisticsAll = aa.outputStatisticsAll;
        
        outputStatistics = [outputStatisticsAll{:}];
        
        if isempty(outputStatistics)
            continue;
        end

        if ~is_GPU
            centrals = outputStatistics(1,:);
            lefts = outputStatistics(2, :);
            rights = outputStatistics(15, :);

            combsLeft = unique([centrals; lefts]', 'rows');
            combsRight = unique([centrals; rights]', 'rows');  %TOO SLOW
        else                
            centrals = gpuArray(outputStatistics(1,:));
            lefts = gpuArray(outputStatistics(2, :));
            rights = gpuArray(outputStatistics(15, :));

            combsLeft = gather(uniqueMyRows(double([centrals; lefts])))';
            combsRight = gather(uniqueMyRows(double([centrals; rights])))';
        end
 
        % convert outputStatistics to the units of the statistical map
        outputStatistics = fromFileToStatMap(outputStatistics, statMapProperties, centreXY, centreZ, centreAngle, offsetConventional, layerID);

        for k = 1:size(combsLeft, 1)
            if isempty(statMapLeft{combsLeft(k,1), combsLeft(k,2)})
                statMapLeft{combsLeft(k,1), combsLeft(k,2)} = uint8(zeros(sizeXY,sizeXY,sizeZ,sizeAngle,sizeAngle,sizeAngle));
            end
        end
        for k = 1:size(combsRight, 1)
            if isempty(statMapRight{combsRight(k,1), combsRight(k,2)})
                statMapRight{combsRight(k,1), combsRight(k,2)} = uint8(zeros(sizeXY,sizeXY,sizeZ,sizeAngle,sizeAngle,sizeAngle));
            end
        end

        for k = 1:size(combsLeft,1)
            tempMatr = uint8(zeros(sizeXY,sizeXY,sizeZ,sizeAngle,sizeAngle,sizeAngle));
            ids = centrals == combsLeft(k,1) & lefts == combsLeft(k,2);
            statisticsTemp = double(outputStatistics(idsLeft, ids));
            inds = sub2ind(size(tempMatr), statisticsTemp(1,:), statisticsTemp(2,:), statisticsTemp(3,:) ,...
                statisticsTemp(4,:), statisticsTemp(5,:), statisticsTemp(6,:));
            t = unique(inds);
            if length(t) > 1
                [freq, idds] = hist(inds, t);
                tempMatr(idds) = freq;
            else
                tempMatr(inds) = length(t);
            end
            statMapLeft{combsLeft(k,1), combsLeft(k,2)} = statMapLeft{combsLeft(k,1), combsLeft(k,2)} + tempMatr;
        end

        for k = 1:size(combsRight,1)
            tempMatr = uint8(zeros(sizeXY,sizeXY,sizeZ,sizeAngle,sizeAngle,sizeAngle));
            ids = centrals == combsRight(k,1) & rights == combsRight(k,2);
            statisticsTemp = double(outputStatistics(idsRight, ids));
            inds = sub2ind(size(tempMatr), statisticsTemp(1,:), statisticsTemp(2,:), statisticsTemp(3,:) ,...
                statisticsTemp(4,:), statisticsTemp(5,:), statisticsTemp(6,:));
            t = unique(inds);
            if length(t) > 1
                [freq, idds] = hist(inds, t);
                tempMatr(idds) = freq;
            else
                tempMatr(inds) = length(t);
            end
            statMapRight{combsRight(k,1), combsRight(k,2)} = statMapRight{combsRight(k,1), combsRight(k,2)} + tempMatr;
        end
    end
end


function combs = uniqueMyRows(matr)
    values = 9999 * matr(1,:) + matr(2,:);
    [~,ia,~] = unique(values);
    combs = matr(:, ia');
end

% unit conversion
function outputStatistics = fromFileToStatMap(outputStatistics, statMapProperties, centreXY, centreZ, centreAngle, offsetConventional, layerID)

    idsAngles = [(6:8)'; (19:21)'];
    
    xyzStep = statMapProperties.xyzStep;
    vectStep = statMapProperties.vectStep;
    
    if layerID == 3 || layerID == 5
        outputStatistics(idsAngles, :) = round(outputStatistics(idsAngles, :)/vectStep) + centreAngle;
        outputStatistics(3, :) = outputStatistics(3, :) - offsetConventional/xyzStep + centreXY;
        outputStatistics(4, :) = outputStatistics(4, :) + centreXY;
        outputStatistics(5, :) = outputStatistics(5, :) + centreZ;
        outputStatistics(16, :) = outputStatistics(16, :) + offsetConventional/xyzStep + centreXY;
        outputStatistics(17, :) = outputStatistics(17, :) + centreXY;
        outputStatistics(18, :) = outputStatistics(18, :) + centreZ;
        
    elseif layerID == 4 || layerID == 6
        outputStatistics(idsAngles, :) = round(outputStatistics(idsAngles, :)/vectStep) + centreAngle;
        outputStatistics(3, :) = outputStatistics(3, :) + centreXY;
        outputStatistics(4, :) = outputStatistics(4, :) - offsetConventional/xyzStep + centreXY;
        outputStatistics(5, :) = outputStatistics(5, :) + centreZ;
        outputStatistics(16, :) = outputStatistics(16, :) + centreXY;
        outputStatistics(17, :) = outputStatistics(17, :) + offsetConventional/xyzStep + centreXY;
        outputStatistics(18, :) = outputStatistics(18, :) + centreZ;
    end
    
end












