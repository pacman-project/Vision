function  [XX, nPrevClusters] = Convert6ToFirstLayer(X, lenCombs, fileForVisualization5Layer)

    load(fileForVisualization5Layer);  %  'triple3OutDepth', 'triple4OutDepth', 'triple5OutDepth'
    nPrevClusters = size(triple5OutDepth, 1);
    
    inds3 = [1,2,6,7,8,9];
    triple3OutDepth = triple3OutDepth(:, inds3);
    
    inds4 = [1,5,6];
    triple4OutDepth = triple4OutDepth(:, inds4);
    triple5OutDepth = triple5OutDepth(:, inds4);  % same inds
    
    
    XX = zeros(lenCombs, 2 * 9 * 9);
    
    for i = 1:lenCombs
        curLine = [];
        
        for mm = 1:3  % triples of the 5th layer
            
            cur5Ind = X(i,mm); % current 5th layer element
            cur5El = triple5OutDepth(cur5Ind, :); 
            
            for k = 1:3 % triples of the 4th layer

                cur4El = triple4OutDepth(cur5El(k), :);  

                for j = 1:3
                    curLine = [curLine, triple3OutDepth(cur4El(j),:)];
                end

            end
        end
        
        XX(i,:) = curLine;
    end

end

