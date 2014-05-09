% this is a function which computes a depth range for each pair of the
% second layer elements

% quant is a quantile which defines depth range (e.g. 0.1)
% output is the following cluster3Depths = zeros(n2Clusters, n2Clusters, lenDisp, 3);
% last parameters are min and max

function [cluster3Depths] = compute3Depths(statistics, n2Clusters, quant, lenDisp)
    
    disp('Analizing depth distribution of pairs ...');
    % initialize the output matrix;
    % the third parameter is displacement
    % the fourth parameter is min(1) max(2) avg(3)
    cluster3Depths = zeros(n2Clusters, n2Clusters, lenDisp, 3); 
    
    % [-127 .. 127] -> [] + 128;
    adder = 128;

    
    firstColumn = statistics(:,1);

    
    for i = 1:n2Clusters
        inds1 = find(firstColumn == i);
        stat = statistics(inds1,:);
        secondColumn = stat(:,2);
        fourthColumn = stat(:,4);
        
        for k = 1:n2Clusters        
            for j = 1:lenDisp             % for each displacement
                if j == 1
                    curColumn = secondColumn;
                    nextCol = 3;
                else
                    curColumn = fourthColumn;
                    nextCol = 5;
                end
                
                inds = find(curColumn == k);
                if isempty(inds)
                    continue;
                end
                
                smallStat = stat(inds, :);
                depths = smallStat(:, nextCol);
                
                depths = double(depths);
                depths(depths < -127) = 127;
                depths(depths > 127)  = 127;
                
                r = length(depths);
                X = zeros(1, 255); 
                for jj = 1:r
                    X(depths(jj)+adder) = X(depths(jj)+adder) + 1;
                end
            
                [depthMin, depthMax, depthAvr] = quantileMy(X, quant, 1-quant, 0.5);

                cluster3Depths(i,k,j,1) = depthMin - adder;
                cluster3Depths(i,k,j,2) = depthMax - adder;
                cluster3Depths(i,k,j,3) = depthAvr - adder;
            end
            
        end
    end
    
end