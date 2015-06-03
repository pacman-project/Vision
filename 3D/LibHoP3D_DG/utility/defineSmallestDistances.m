function [smallestDist, smIDs] = defineSmallestDistances(distances, numEls)

    [distances, idx] = sort(distances, 'ascend');
    
    while distances(numEls) == distances(numEls+1)
        numEls = numEls + 1;
    end
    
    smIDs = idx(1:numEls);
    smallestDist = distances(1:numEls);
end

