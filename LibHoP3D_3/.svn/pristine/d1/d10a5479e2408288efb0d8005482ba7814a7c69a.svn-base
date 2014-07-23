% this is the function to present a line in terms of parts
% In contrast to the previous version, this function rely on optimization

% we optimize:  min [error + wOverlap * overlap - wCoverage*coverage ]
% error = 

% wCoverage - weight parameter in the optimization function (above),
% wOverlap - weight of the overlap term in the above optimization function
% errors - errrs of each element (distances to the nearest vocabulary element)
% overlaps - overlaps between candidate elements. 
% coverage - coverage of each candidate element.

function [w] = LineDiscretizationOpt(wCoverage, wOverlap, errors, overlaps, coverages, combs)
  
    len = length(errors);
    
    errMatrix = zeros(len, len);
    covMatrix = zeros(len, len);
    overMatrix = zeros(len, len);
        
    % fill these matrices
    for i = 1:len
        errMatrix(i,i) = errors(i);
        covMatrix(i,i) = coverages(i);
        if i == 1
            errMatrix(1,2) = -errors(i)/4;
            covMatrix(1,2) = -coverages(i)/4;
            overMatrix(1,2)= overlaps(i);
        elseif i == len
            errMatrix(i, i-1) = -errors(i)/4;
            covMatrix(i, i-1) = -coverages(i)/4;
            overMatrix(i, i-1)= overlaps(i-1);
        else % all the other raws
            errMatrix(i,i+1) = -errors(i)/4;
            errMatrix(i,i-1) = -errors(i)/4;
            covMatrix(i,i+1) = -coverages(i)/4;
            covMatrix(i,i-1) = -coverages(i)/4;
            overMatrix(i,i+1) = overlaps(i);
            overMatrix(i,i-1) = overlaps(i-1);
        end
    end
    
    minOpt = 10^6;
    w = zeros(1, len); % weight vector
    
    if isempty(combs) % generate combinations right here
        
        numComb = 2^len; % number of possible combinations
        combStart = 2^(len-2);

        for i = combStart:numComb
             str = dec2bin(i-1, len);
             % check whether this is a good combination
             k = strfind(str, '000'); % features for bad combinations
             if isempty(k)
                 kk = strfind(str, '11111');
                 if isempty(kk) % combination looks good       
                     curComb = zeros(1,len);
                     for j = 1:len
                         curComb(j) = str2double(str(j));
                     end

                     curW = curComb;
                     % vector curW is now ready. Compute optimization function
                     % [error + wCoverage*coverage - wOverlap * overlap]
                     curOpt = curW*(errMatrix +  wOverlap*overMatrix - wCoverage*covMatrix)*curW';
                     if curOpt < minOpt
                         minOpt = curOpt;
                         w = curW;
                     end
                 end
             end
        end
    
    else
        [r,c] = size(combs);
        for i = 1:r
             curW = combs(i,:);
             % vector curW is now ready. Compute optimization function
             % [error + wCoverage*coverage - wOverlap * overlap]
             curOpt = curW*(errMatrix +  wOverlap*overMatrix - wCoverage*covMatrix)*curW';
             if curOpt < minOpt
                 minOpt = curOpt;
                 w = curW;
             end
        end
        
    end
end