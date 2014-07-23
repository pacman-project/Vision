% this function is used to replace data as follows
    %      9 5 8      
    %      2 1 4      
    %      6 3 7  --->   [central, top, top_depth, bottom, bottom_depth]
    
% is_exact shows whether we want an exact translation

function [statistics, curTS] = replaceTripleWith3Elements(outputStatistics, curTS, triples3Out, n2Clusters, nClusters, is_exact)
    
    a = [9,5,8; 2,1,4; 6,3,7];
    len = 0;
    statistics = zeros(curTS, 5);
    inds = [];
    
    if is_exact
        % create a table for faster matching
        table =  zeros(n2Clusters, n2Clusters, n2Clusters);
        n3Clusters = size(triples3Out ,1);
        for k = 1:n3Clusters
            cur = triples3Out(k,:);
            ind = zeros(3,1);
            ind(1) = compute2elementIndex(cur(1), cur(2), nClusters);
            ind(2) = compute2elementIndex(cur(3), cur(4), nClusters);
            ind(3) = compute2elementIndex(cur(5), cur(6), nClusters);
            table(ind(1), ind(2), ind(3)) = k; 
        end

        for i = 1:curTS
            els3 = zeros(1,3);
            is_ok = true;
            curLine = outputStatistics(i,:);
            for j = 1:3
                cur = a(j,:);   % ex.  9,5,8
                if j ~= 2
                    cols = [(cur(1)-1)*2, (cur(2)-1)*2, (cur(3)-1)*2]; % ex.  9,5,8
                else            % ex.   2,1,4
                    cols = [(cur(1)-1)*2, 1, (cur(3)-1)*2]; % ex.  9,5,8
                end
                curStat = curLine(cols);
                % check if it belongs to the dictionary
                num = table(curStat(1), curStat(2), curStat(3));
                els3(j) = num;
                if num == 0 % element is not in the vocabulary
                    is_ok = false;
                    break;
                end
            end
            
            if is_ok
                % compute relative depths of the top and bottom elements
                depthTop = curLine(3);
                depthBottom = curLine(7);
                len = len + 1;
                statistics(len,:) = [els3(2), els3(1), depthTop, els3(3), depthBottom];
            end

        end
        statistics = statistics(1:len,:);
        curTS = len;
    else
        disp('ERROR. This is not implemented yet');
    end
    

end



