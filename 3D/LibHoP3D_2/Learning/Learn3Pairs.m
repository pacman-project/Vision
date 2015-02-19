% this is a function for learning pairs of the third layer based on statistics
% statistics has the following format:
% central, left, depthLeft, right, depthRight

% Output:
% tablePairs = zeros(2, nPrevClusters, nPrevClusters)
% 2 - left right
% nPrevClusters - central element
% nPrevClusters - another element


function [tablePairs, numPairs] = Learn3Pairs(statistics, curTS, nPrevClusters, thresh3Pair, lenDisp)

    disp('Learning pairs of the third layer ...');
    tablePairs = zeros(lenDisp, nPrevClusters, nPrevClusters);
    
    
    parfor i = 1:lenDisp % for each displacement
        
        tableTemp = zeros(nPrevClusters, nPrevClusters);
        for j = 1:curTS
            tableTemp(statistics(j,1), statistics(j,2*i)) = tableTemp(statistics(j,1), statistics(j,2*i)) + 1;  % increase the counter
        end
        
        tableTemp(tableTemp < thresh3Pair) = 0;
        tablePairs(i, :, :) = tableTemp;
    end

    numPairs = length(tablePairs(tablePairs > 0));
    
end