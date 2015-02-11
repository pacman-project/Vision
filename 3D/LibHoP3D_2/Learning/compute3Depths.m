% this is a function which computes a depth range for each pair of the
% second layer elements

% quant is a quantile which defines depth range (e.g. 0.1)
% output is the following cluster3Depths = zeros(nPrevClusters, nPrevClusters, numDisps, 3);
% last parameters are min and max

function [cluster3Depths] = compute3Depths(statistics, nPrevClusters, quant, numDisps, maxRelDepth, is_GPU_USED)
 
    
    disp('Analizing depth distribution of pairs ...');
    % initialize the output matrix;
    % the third parameter is displacement
    % the fourth parameter is min(1) max(2) avg(3)
    cluster3Depths = zeros(nPrevClusters, nPrevClusters, numDisps, 3);
    
    %          [-127 .. 127]       -> [] + 128;
    % becomes: [- dThresh, dThresh] = [] + addder;
    adder = maxRelDepth + 1;
    dThresh = maxRelDepth;
    
    if is_GPU_USED
        statistics = gpuArray(statistics);
    end
    
    histEdges = 1 : (2 * dThresh + 1);
    
    for i = 1:nPrevClusters
        
        inds = statistics(:,1) == i;
        if max(inds) == 0
            continue;
        end

        stat = statistics(inds, :);
        
        for j = 1:numDisps 
            curColumn = j*2;
        
            for k = 1:nPrevClusters
                indSS = stat(:,curColumn) == k;
                
                if max(indSS) ~= 0
                    smallStat = stat(indSS, :);
                    % define list of relative depths
                    if is_GPU_USED
                        depths = gather(smallStat(:, curColumn + 1));
                    else
                        depths = smallStat(:, curColumn + 1);
                    end
                    
                    depths = double(depths);
                    depths(depths < -dThresh) = -dThresh;
                    depths(depths > dThresh)  =  dThresh;  
                    depths = depths + adder;
                    X = histc(depths, histEdges);
                    
%                     for jj = 1:r
%                         X(depths(jj)+adder) = X(depths(jj)+adder) + 1;
%                     end

                    [depthMin, depthMax, depthAvr] = quantileMy(X, quant, 1-quant, 0.5);

                    cluster3Depths(i,k,j,1) = depthMin - adder;
                    cluster3Depths(i,k,j,2) = depthMax - adder;
                    cluster3Depths(i,k,j,3) = depthAvr - adder;
                end

            end
        end
    end

    
end




%     if (nargin == 4)
%         maxRelDepth = 127;
%     end 
%     
%     disp('Analizing depth distribution of pairs ...');
%     % initialize the output matrix;
%     % the third parameter is displacement
%     % the fourth parameter is min(1) max(2) avg(3)
%     cluster3Depths = zeros(nPrevClusters, nPrevClusters, numDisps, 3); 
%     
%     %          [-127 .. 127]       -> [] + 128;
%     % becomes: [- dThresh, dThresh] = [] + addder;
%     adder = maxRelDepth + 1;
%     dThresh = maxRelDepth;
% 
%     
%     firstColumn = statistics(:,1);
% 
%     
%     for i = 1:nPrevClusters
%         inds1 = firstColumn == i;
%         stat = statistics(inds1,:);
%         secondColumn = stat(:,2);
%         fourthColumn = stat(:,4);
%         
%         for k = 1:nPrevClusters        
%             for j = 1:numDisps             % for each displacement
%                 if j == 1
%                     curColumn = secondColumn;
%                     nextCol = 3;
%                 elseif j == 2
%                     curColumn = fourthColumn;
%                     nextCol = 5;
%                 end
%                 
%                 inds = find(curColumn == k);
%                 if isempty(inds)
%                     continue;
%                 end
%                 
%                 smallStat = stat(inds, :);
%                 depths = smallStat(:, nextCol);
%                 
%                 depths = double(depths);
%                 depths(depths < -dThresh) = -dThresh;
%                 depths(depths > dThresh)  =  dThresh;
%                 
%                 r = length(depths);
%                 X = zeros(1, 2 * dThresh + 1); 
%                 for jj = 1:r
%                     X(depths(jj)+adder) = X(depths(jj)+adder) + 1;
%                 end
%             
%                 [depthMin, depthMax, depthAvr] = quantileMy(X, quant, 1-quant, 0.5);
% 
%                 cluster3Depths(i,k,j,1) = depthMin - adder;
%                 cluster3Depths(i,k,j,2) = depthMax - adder;
%                 cluster3Depths(i,k,j,3) = depthAvr - adder;
%             end
%             
%         end
%     end
%     
% end