% this is to compute score for the optimization function

% we penalize distances less than thresh

function [score] = computeScore(X, frequencies, adopted, alpha, betta, distThresh)
    
    [rX, ~] = size(X);
    lenAdopted = length(adopted);
    
    % distances is n*rX matrix
    distances = Isodata_distances(X, adopted, rX, lenAdopted, false, false);
    
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

