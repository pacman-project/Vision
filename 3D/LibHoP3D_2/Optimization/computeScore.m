% this is to compute score for the optimization function

% we penalize distances less than thresh

function [score] = computeScore(X, frequencies, adopted, added, alpha, betta, distThresh, distXtoInitAdopted)
    
    rX = size(X, 1);  
    lenAdded = size(added, 1);
    distances = distXtoInitAdopted;
    
    % distances is rX*lenAdded  matrix
    if lenAdded > 0
        distancesAdd = Isodata_distances(X, added, rX, lenAdded, false, false);
        distances = [distances, distancesAdd];
    end
    
    lenAdopted = length(adopted) + lenAdded;
    
    score = betta * lenAdopted;
    minDists = min(distances, [], 2);
    
    for i = 1:rX
        score = score + frequencies(i) * minDists(i);
    end
    
    if alpha ~= 0
        withinDist = Isodata_distances(adopted, adopted, lenAdopted, lenAdopted, false, false); 
        score = score + sum(sum(withinDist<=distThresh)) * alpha;
        score = score + sum(sum(withinDist<=distThresh/2)) * alpha;  % small distances are penalized twice
    end
    
end

