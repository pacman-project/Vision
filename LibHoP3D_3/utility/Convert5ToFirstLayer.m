function  [XX, nPrevClusters] = Convert5ToFirstLayer(X, lenCombs, fileForVisualization4Layer)

    load(fileForVisualization4Layer);  %  'triple3OutDepth', 'triple4OutDepth'
     nPrevClusters = size(triple4OutDepth, 1);
    
    inds3 = [1,2,6,7,8,9];
    triple3OutDepth = triple3OutDepth(:, inds3);
    
    inds4 = [1,5,6];
    triple4OutDepth = triple4OutDepth(:, inds4);
    
    
    XX = zeros(lenCombs, 18 * 3);
    
    for i = 1:lenCombs
        curLine = [];
        
        for k = 1:3
            
            cur4Ind = X(i,k); % current 4th layer element
            cur4El = triple4OutDepth(cur4Ind, :);  % in format: top central bottom
            
            for j = 1:3
                curLine = [curLine, triple3OutDepth(cur4El(j),:)];
            end
            
        end
        
        XX(i,:) = curLine;
    end

end

