% this is the script for greedy part selection for the layer 3.

function [partsOut, coverageOut, lenOut] = partSelectionNew(nClusters, n2Clusters, statistics3LayerSieved, statistics3LayerAggregated, ...
                             dataSetNumber, fieldSize, list_depth, lenF, ...
                             iterations, layerID, fileForVisualizationPrevLayer, displ3, displ5, displ7, cluster1Centres, ...
                             depthStep, lenSelected, numSimilar, is_GPU_USED)
                                                 
    % output variables
    coverageOut = zeros(iterations, 1);
    partsOut = zeros(iterations,3);
    lenOut = 0;
                           
    halfFieldSize = floor(fieldSize/2);    %         for example fieldSize = [17, 5, 71];
    fieldCenter = halfFieldSize + 1;
    
    dx = halfFieldSize(1);
    dy = halfFieldSize(2);
       
    a = load(statistics3LayerSieved);          %       'statistics', 'clusterCurDepths', 'outputCoords'
    b = load(statistics3LayerAggregated);      %       'X' ,'frequencies', 'curTS', 'triples'
    
    statistics = a.statistics;
    
    lenStat =  size(statistics, 1);
    clusterCurDepths = a.clusterCurDepths;
    outputCoords = a.outputCoords;
    X = b.X;
    curTS = b.curTS;
    
    clear('a', 'b');
    
    lenCombs = size(X, 1);
   
    downsamplingScheme = 3;
    
    [X_first, nPrevClusters, emptyIndicator] = convertToSurfaceDescriptor(X, lenCombs, layerID, nClusters, n2Clusters, fileForVisualizationPrevLayer{layerID-1},  ...
                       displ3, displ5, displ7, fieldCenter, cluster1Centres, downsamplingScheme, clusterCurDepths);  % this is a first layer descriptor
                   
    X_first(emptyIndicator == 1) = fieldSize(3)*3*depthStep;
    coverage = zeros(lenCombs, 1);      % just something to initialize
    coverage = double(coverage);
    
    % this is done to speed up a function of reconstruction
    
    if nPrevClusters > 800
        is_sparse = true;
    else
        is_sparse = false;
    end
    
    % create a matrix once
    if is_sparse
        parfor i = 1:nPrevClusters
            triples{i} = sparse(zeros(nPrevClusters + 1, nPrevClusters + 1));
        end
    else
        triples = [];  % otherwise create a matrix in the function recomputeCoverage
    end 
    
    indsLine = [2,1,4];
    stat = statistics(:, indsLine);
    
    [coverage] = recomputeCoverage(statistics, outputCoords, coverage, X, 1, dx, dy,...
                                    list_depth, lenF, nPrevClusters + 1, lenCombs, is_sparse, triples);
                                
%     coverage(1:kk-1) = ones(1:(kk-1), 1)*(max(coverage)+10);
%     coverageOut(1:kk-1) = ones(1:(kk-1), 1)*(max(coverage)+10);

 
    % sort parts according to their coverage
    [coverage, inds] = sort(coverage, 'descend');
    X = X(inds, :);
    X_first = X_first(inds, :);   
    X = int16(X);
    
    if is_GPU_USED
        stat = gpuArray(stat);
    end
    
    for i = 1:iterations   %kk:iterations+kk   % for each part 
        
        str = ['element - ', num2str(i)];
        disp(str);
        str = ['coverage - ', num2str(coverage(i))];
        disp(str);
        
        if coverage(i) == 0
            return;  % all parts are selected
        end 
        
        % merge the selected elements with all elements with distance
        % less than meargeThresh
        
        distance = Integral_distances(X_first, X_first(i,:), lenCombs, 1, false, false);
        [smallestDists, smIDs] = defineSmallestDistances(distance, numSimilar);
        
        
        curPart = X(smIDs, :); 
        numSimilarUpd = size(curPart, 1); % may be a bit different from the previous one
        
        % find position of curParts in all images
        if is_GPU_USED

            inds = arrayfun(@isMemberMy, stat(:,1), stat(:,2), stat(:,3));
            elPositions = outputCoords(inds, :);
        else
            
            inds = ismember(stat, curPart, 'rows');
            elPositions = outputCoords(inds, :);
        end
      
        % project this part to the images
        % curElPosition is already sorted by image number
        firstColumn = elPositions(:,1);
        
        parfor j = 1:lenF
            
            % extract all related to the image
            % 1) find the current position in the elPositions
            ind = find(firstColumn == j);
            lenE = length(ind);
            
            if lenE > 0
                % open the image
                I = imread(list_depth{j});
                [r,c, ch] = size(I);

                % the right image is now open
                for k = 1:lenE  % fill all positions for this image

                    x = elPositions(ind(k), 2);
                    y = elPositions(ind(k), 3);
                    
                    indStartX = x-dx;
                    indEndX = x+dx;
                    indStartY = y-dy;
                    indEndY = y+dy;
                    
                    if indStartX < 1
                        indStartX = 1;
                    elseif indEndX > c
                        indEndX = c;
                    end
                    
                    if indStartY < 1
                        indStartY = 1;
                    elseif indEndY > r
                        indEndY = r;
                    end        
                    
                    I(indStartY:indEndY ,indStartX:indEndX, 3) = zeros(2*dy+1, 2*dx+1); 
                end

                % save the last image
                I = uint16(I);
                imwrite(I, list_depth{j}, 'png');
            end
            
        end
        
        coverageOut(i) = coverage(i);
        partsOut(i, :) = X(i,:);
        lenOut = lenOut + 1;
        
 
        %------------------------------------------------------------------
        % update coverage of the following parts
        lenRec = 30;
        
        if mod(i, 20) == 0 || (coverage(i+1)<= lenSelected) % full recompute
            lenRec = lenCombs - i;
        elseif mod(i, 8) == 0
            lenRec = 200;
        elseif mod(i, 3) == 0
            lenRec = 60;
        end

        [coverage] = recomputeCoverage(statistics, outputCoords, coverage, X, i+1, dx, dy,...
                                        list_depth, lenF, nPrevClusters + 1, lenRec, is_sparse, triples);
        

        % sort parts according to their coverage (sort only the rest of the array)
        [coverage(i+1:end), inds] = sort(coverage(i+1:end), 'descend');
        inds = inds+i;
        X(i+1:end,:) = X(inds, :);
        X_first(i+1:end,:) = X_first(inds, :);
        
        if (mod(i, 10) == 0) && (coverage(i+1)<= lenSelected)  % maximal remaining coverage is less than we expect
            return;
        end
            
        save('Temp/partsOut.mat', 'partsOut');
    end
    
    
    % nested function impemented on gpu
    function isTrue = isMemberMy(col1, col2, col3)
        isTrue = false;
        for jj = 1:numSimilarUpd
            isTrue = isTrue | (curPart(jj,1) == col1 & curPart(jj,2) == col2 & curPart(jj,3) == col3);
        end
        
    end

    if is_GPU_USED
        reset(gpuDevice(1));
    end
       
end
    





