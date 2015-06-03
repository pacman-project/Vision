% this is the function to present a line in terms of parts
% In contrast to the previous version, this function rely on optimization

% we optimize:  min [error + wOverlap * overlap - wCoverage*coverage ]
% error = 

% wCoverage - weight parameter in the optimization function (above),
% wOverlap - weight of the overlap term in the above optimization function
% errors - errrs of each element (distances to the nearest vocabulary element)
% overlaps - overlaps between candidate elements. 
% coverage - coverage of each candidate element.

function [w] = LineInhibitionOpt(wCoverage, wOverlap, errors, overlapsClose, overlapsFar, coverages, combs)
  
    len = length(errors);
    if len == 1
        w = 1;
        return;
    elseif len == 2
        w = [1 0];
        return;
    end
    
    errMatrix = zeros(len, len);
    covMatrix = zeros(len, len);
    overMatrix = zeros(len, len);
    
    errorMultiplier1 = 0.3;
    errorMultiplier2 = 0.1;
    
    coverageSubtracter1 = 1.5;
    coverageSubtracter2 = 0.5;
    
        
    % fill these matrices
    for i = 1:len
        errMatrix(i,i) = errors(i);
        covMatrix(i,i) = coverages(i);
        if i == 1
            errMatrix(1,2) = -errors(i)*errorMultiplier1;
            errMatrix(1,3) = -errors(i)*errorMultiplier2;
            if coverages(i) ~= 0
                covMatrix(1,2) = -coverageSubtracter1;
                covMatrix(1,3) = -coverageSubtracter2;
            end
            overMatrix(1,2)= overlapsClose(i);
            overMatrix(1,3)= overlapsFar(i);
        elseif i == len
            errMatrix(i, i-1) = -errors(i)*errorMultiplier1;
            errMatrix(i, i-2) = -errors(i)*errorMultiplier2;
            if coverages(i) ~= 0
                covMatrix(i, i-1) = -coverageSubtracter1;
                covMatrix(i, i-2) = -coverageSubtracter2;
            end
            overMatrix(i, i-1)= overlapsClose(i-1);
            overMatrix(i, i-2)= overlapsFar(i-2);
        else % all the other raws
            errMatrix(i,i+1) = -errors(i)*errorMultiplier1;
            errMatrix(i,i-1) = -errors(i)*errorMultiplier1;
            if coverages(i) ~= 0
                covMatrix(i,i+1) = -coverageSubtracter1;
                covMatrix(i,i-1) = -coverageSubtracter1;
            end
            overMatrix(i,i+1) = overlapsClose(i);
            overMatrix(i,i-1) = overlapsClose(i-1);
            
            if i ~= 2
                overMatrix(i,i-2) = overlapsFar(i-2);
                errMatrix(i, i-2) = -errors(i)*errorMultiplier2;
                if coverages(i) ~= 0
                    covMatrix(i, i-2) = -coverageSubtracter2;
                end
            end
            if i ~= len - 1
                overMatrix(i,i+2) = overlapsFar(i);
                errMatrix(i, i+2) = -errors(i)*errorMultiplier2;
                if coverages(i) ~= 0
                    covMatrix(i, i+2) = -coverageSubtracter2;
                end
            end
        end
    end
    
    minOpt = 10^6;
    w = zeros(1, len); % weight vector
    
    if isempty(combs) % generate combinations right here
        
        % call a function to generate combinations
        [combs] = createCombinations(len);
    end
    
    % apply combinations to find a minimizer of a quadratic form
    
    [r,~] = size(combs);
    for i = 1:r
         % vector curW is now ready. Compute optimization function
         % [error + wCoverage*coverage - wOverlap * overlap]
         curOpt = combs(i,:)*(errMatrix +  wOverlap*overMatrix - wCoverage*covMatrix)*combs(i,:)';
         if curOpt < minOpt
             minOpt = curOpt;
             w = combs(i,:);
         end
    end

    
    a = w;
end