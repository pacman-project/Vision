function  [XX, nPrevClusters] = Convert4ToFirstLayer(X, lenCombs, fileForVisualization3Layer)

    load(fileForVisualization3Layer);  % triple3OutDepth
    nPrevClusters = size(triple3OutDepth, 1);
    
    inds = [1,2,6,7,8,9];
    
    triple3OutDepth = triple3OutDepth(:, inds);
    XX = zeros(lenCombs, 18);
    
    for i = 1:lenCombs
        curLine = [];
        
        for j = 1:3
            curLine = [curLine, triple3OutDepth(X(i,j),:)];
        end
        
        XX(i,:) = curLine;
    end

end

