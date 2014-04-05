% this is for the line discretization of the first layer
% this function is based on LineDiscretizationOpt

% fx - line of filter responces
% nClusters, thresh, clusterSize, clusterCentres - parameters of the first layer

% nearestClusters == 0 means there is no data in these locations

function [outLine, outErrors] = lineDiscretizationOptLayer1(nearestClusters, errors, wCoverage, wOverlap, combs, largestLine)
    
    strLen = length(nearestClusters);
    coverages = ones(1, strLen) * 6;
    coverages(nearestClusters == 0) = 0;
    maxOverlap = 2;
    overlaps = zeros(1, strLen);
    
    for i = 1:strLen-1
        if nearestClusters(i) == nearestClusters(i+1)
            overlaps(i) = maxOverlap;
        else
            overlaps(i) = maxOverlap/4;  % we can do something even better
        end
    end
    overlaps(strLen) = 2; % seems to be not used
      
    % decompose if line is too long
    els = [];
    ers = [];
    if strLen > largestLine % do decomposition
        while strLen >= 2 * largestLine % feed it by portions of size largestLine
            curEls = LineDiscretizationOpt(wCoverage, wOverlap, errors(1:largestLine), overlaps(1:largestLine-1), coverages(1:largestLine), combs);
            curErs = errors(1:largestLine) .* curEls; 
            els = [els, curEls];
            ers = [ers, curErs];
            errors = errors(largestLine + 1:end);
            overlaps = overlaps(largestLine + 1:end); % should be ok
            coverages = coverages(largestLine + 1:end);            
            strLen = length(errors);
            if els(end) == 1
                coverages(1) = coverages(1) - 2;
            end
        end
        
        if (strLen < 2 * largestLine)  && (strLen > largestLine)
            tempLen = ceil(strLen/2);
            curEls = LineDiscretizationOpt(wCoverage * 1.3, wOverlap, errors(1:tempLen), overlaps(1:tempLen-1), coverages(1:tempLen), []);
            curErs = errors(1:tempLen) .* curEls;
            els = [els, curEls];
            ers = [ers, curErs];
            errors = errors(tempLen + 1:end);
            overlaps = overlaps(tempLen + 1:end); % should be ok
            coverages = coverages(tempLen + 1:end);
            if els(end) == 1
                coverages(1) = coverages(1) - 2;
            end
        end
    end

    [curEls] = LineDiscretizationOpt(wCoverage*1.5, wOverlap, errors, overlaps(1:end-1), coverages, []); % feed the rest
    curErs = errors .* curEls;
    els = [els, curEls]; % this is 0-1 line  [0,1,1,0,1,...]
    ers = [ers, curErs];

    outLine = els.* nearestClusters;  % make it [0, 4, 5, 0, 2,...]
    ers(nearestClusters == 0) = 0;
    outErrors = ers;
end















