% this is a function for learning pairs of the third layer based on statistics
% statistics has the following format:
% central, left, depthLeft, right, depthRight

% Output:
% tablePairs = zeros(2, n2Clusters, n2Clusters)
% 2 - left right
% n2clusters - central element
% n2clusters - another element


function [tablePairs, numPairs] = Learn3Pairs(statistics, curTS, n2Clusters, thresh3Pair, lenDisp)

    disp('Learning pairs of the third layer ...');
    tablePairs = zeros(lenDisp, n2Clusters, n2Clusters);
    
    for i = 1:lenDisp % for each displacement
        % extract subset
        indC = [1, 2*i, 2*i+1];
        curStat = statistics(:, indC);
        
        % 1) measure frequency of each pair

        for j = 1:curTS
            tablePairs(i, curStat(j,1), curStat(j,2)) = tablePairs(i, curStat(j,1), curStat(j,2)) + 1;  % increase the counter
        end
        
        % 2) filter out the least frequent pairs
        for jj = 1:n2Clusters
            for ii = 1:n2Clusters
                if tablePairs(i, ii, jj) < thresh3Pair
                    tablePairs(i, ii, jj) = 0; 
                end        
            end
        end
        
        i
    end
    
    numPairs = length(tablePairs(tablePairs > 0));
end