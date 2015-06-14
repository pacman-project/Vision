% Computes a table of elements to be replaced according to pre-defined
% grammar rules

% abstractionLevel == 0   no abstraction
% abstractionLevel == 1   replace convex parts with less convex
% abstractionLevel == 2   replace convex parts with less convex And twisted
% parts with less twisted

function [abstractionTable] = computeAbstractionTable(X, nClusters, n2Clusters, abstractionLevel, table2)


    disp('compute the abstraction table...');
    lenCombs = size(X, 1);
    abstractionTable = zeros(lenCombs, 1);
    
    triples = zeros(n2Clusters,n2Clusters,n2Clusters);
    
    % fill a table triples
    
    for i = 1:lenCombs
        triples(X(i, 1), X(i, 2), X(i, 3)) = i;
    end
    
    
    if abstractionLevel == 1 || abstractionLevel == 2
        for i = i:lenCombs
            curEl = X(i,:);

            left   = curEl(1);
            centre = curEl(2);
            right  = curEl(3);
            changed = false;

            if left == right && left ~= centre % something wiered in the middle of the part
                if abs(centre - right) == 1 || abs(centre - right) == nClusters  % not too wiered
                    centre = right;

                    indNew = triples(left, centre, right);
                    abstractionTable(i) = indNew;            
                    continue;
                end
            end

            if left ~= centre || right~=centre % non flat element (abstraction may be applied)

                clusterXL = table2(left, 1);
                clusterYL = table2(left, 2);
                clusterXR = table2(right, 1);
                clusterYR = table2(right, 2);
                clusterXC = table2(centre, 1);
                clusterYC = table2(centre, 2); 

                if clusterXL == clusterXC-2  || clusterXL == clusterXC-3
                    clusterXL = clusterXC-1;
                    changed = true;
                elseif clusterXL == clusterXC+2 || clusterXL == clusterXC+3
                    clusterXL = clusterXC+1;
                    changed = true;
                end

                if clusterXR == clusterXC-2 || clusterXR == clusterXC-3
                    clusterXR = clusterXC-1;
                    changed = true;
                elseif clusterXR == clusterXC+2 || clusterXR == clusterXC+3
                    clusterXR = clusterXC+1;
                    changed = true;
                end

                if abstractionLevel == 2 % replace twisted parts with less twisted ones

                    if clusterYL == clusterYC-2 || clusterYL == clusterYC-3
                        clusterYL = clusterYC-1;
                        changed = true;
                    elseif clusterYL == clusterYC+2 || clusterYL == clusterYC+3
                        clusterYL = clusterYC+1;
                        changed = true;
                    end

                    if clusterYR == clusterYC-2 || clusterYR == clusterYC-3
                        clusterYR = clusterYC-1;
                        changed = true;
                    elseif clusterYR == clusterYC+2 || clusterYR == clusterYC+3
                        clusterYR = clusterYC+1;
                        changed = true;
                    end

                end

                % Recover the index of new element
                if changed
                    leftNew = compute2elementIndex(clusterXL, clusterYL, nClusters);
                    rightNew = compute2elementIndex(clusterXR, clusterYR, nClusters);

                    % find the new element in the table triples
                    indNew = triples(leftNew, centre, rightNew);
                    abstractionTable(i) = indNew;
                end
            end
        end

    end
end   




