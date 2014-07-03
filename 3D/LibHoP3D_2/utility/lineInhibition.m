% this is for inhibition of the first layer elements
% this function is based on lineInhibitionOpt

% fx - line of filter responces
% nClusters, thresh, clusterSize, clusterCentres - parameters of the first layer

% nearestClusters == 0 means there is no data in these locations

function [outLine, outErrors] = lineInhibition(nearestClusters, errors, wCoverage, wOverlap, combs, largestLine)
    
    elementSize = 5;

    strLen = length(nearestClusters);
    coverages = ones(1, strLen) * elementSize;
    coverages(1) = elementSize - 2;
    coverages(end) = elementSize - 2;
    coverages(nearestClusters == 0) = 0;
    
    maxOverlap = 3;
    minOverlap = 1;
    
    scaleFactor = 2;
    
    overlapsClose = zeros(1, strLen);
    overlapsFar = zeros(1, strLen);
    
    %% This is done to reduce overlap of the closest elements
    
    for i = 1:strLen-1
        if nearestClusters(i) ~= 0 && nearestClusters(i+1) ~= 0
            if nearestClusters(i) == nearestClusters(i+1) 
                overlapsClose(i) = maxOverlap;
            else
                overlapsClose(i) = maxOverlap/scaleFactor;
            end
        end
    end
    
    for i = 1:strLen-2
        if nearestClusters(i) ~= 0 && nearestClusters(i+2) ~= 0
            if nearestClusters(i) == nearestClusters(i+2)
                overlapsFar(i) = minOverlap;
            else
                overlapsFar(i) = minOverlap/scaleFactor;
            end
        end
    end
    
%     overlapsClose(strLen) = maxOverlap; % last element
%     overlapsFar(strLen-1) = minOverlap; % last element
%     overlapsFar(strLen)   = minOverlap; % last element  % NOT USED 
    
    %%
      
    % decompose if line is too long
    els = [];
    ers = [];
    if strLen > largestLine % do decomposition
        while strLen >= 2 * largestLine % feed it by portions of size largestLine
            curEls = LineInhibitionOpt(wCoverage, wOverlap, errors(1:largestLine), overlapsClose(1:largestLine-1), overlapsFar(1:largestLine-1), coverages(1:largestLine), combs);
            curErs = errors(1:largestLine) .* curEls; 
            els = [els, curEls];
            ers = [ers, curErs];
            errors = errors(largestLine + 1:end);
            overlapsClose = overlapsClose(largestLine + 1:end); % should be ok
            overlapsFar = overlapsFar(largestLine + 1:end); % should be ok
            coverages = coverages(largestLine + 1:end);            
            strLen = length(errors);
            if els(end) == 1
                coverages(1) = coverages(1) - maxOverlap;
                coverages(2) = coverages(2) - minOverlap;
                overlapsClose(1) = overlapsClose(1) + 4;
                overlapsFar(1) = overlapsFar(1) + 2;
                errors(1) = 100;  % not good
                errors(2) = errors(2)*2;  % not good
            end
        end
        
        if (strLen < 2 * largestLine)  && (strLen > largestLine)
            tempLen = ceil(strLen/2);
            curEls = LineInhibitionOpt(wCoverage * 1.3, wOverlap, errors(1:tempLen), overlapsClose(1:tempLen-1), overlapsFar(1:tempLen-1), coverages(1:tempLen), []);
            curErs = errors(1:tempLen) .* curEls;
            els = [els, curEls];
            ers = [ers, curErs];
            errors = errors(tempLen + 1:end);
            overlapsClose = overlapsClose(tempLen + 1:end); % should be ok
            overlapsFar = overlapsFar(tempLen + 1:end); % should be ok
            coverages = coverages(tempLen + 1:end);
            if els(end) == 1
                coverages(1) = coverages(1) - maxOverlap;
                coverages(2) = coverages(2) - minOverlap;
                overlapsClose(1) = overlapsClose(1) + 4;
                overlapsFar(1) = overlapsFar(1) + 2;
                errors(1) = 100;  % not good
                errors(2) = errors(2)*2;  % not good
            end
        end
    end

    [curEls] = LineInhibitionOpt(wCoverage * 1.5, wOverlap, errors, overlapsClose(1:end-1), overlapsFar(1:end-1), coverages, []); % feed the rest
    curErs = errors .* curEls;
    els = [els, curEls]; % this is 0-1 line  [0,1,1,0,1,...]
    ers = [ers, curErs];

    outLine = els.* nearestClusters;  % make it [0, 4, 5, 0, 2,...]
    ers(nearestClusters == 0) = 0;
    outErrors = ers;
end















