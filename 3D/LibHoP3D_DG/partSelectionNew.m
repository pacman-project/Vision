% this is the script for greedy part selection for the layer 3.

function [partsOut, coverageOut, lenOut] = partSelectionNew(nClusters, n2Clusters, statistics3LayerSieved, statistics3LayerAggregated, ...
                             dataSetNumber, fieldSize, list_depth, lenF, iterations, layerID, fileForVisualizationPrevLayer, ...
                             cluster1Centres, lenSelected, numSimilar, offsetsConventional, depthStep)
                                                 
    % output variables
    coverageOut = zeros(iterations, 1);
    partsOut = zeros(iterations,3);
    lenOut = 0;
                           
    halfFieldSize = ceil(fieldSize/2);    %         for example fieldSize = [17, 5, 71];
    fieldCenter = halfFieldSize;
    
    dx = halfFieldSize(1);
    dy = halfFieldSize(2);
    
    
    a = load(statistics3LayerSieved);          %       'statistics', 'clusterCurDepths', 'outputCoords'
    b = load(statistics3LayerAggregated);      %       'X' ,'frequencies', 'curTS', 'triples'
    
    statistics = a.statistics;
    
    lenStat =  size(statistics, 1);
    clusterCurDepths = a.clusterCurDepths;
    outputCoords = a.outputCoords;
    outputScales = a.outputScales;
    outputFrames = a.outputFrames;
     
    X = b.X;
    curTS = b.curTS;
    
    clear('a', 'b');
    
    lenCombs = size(X, 1);
   
    downsamplingScheme = 3;
    
    [X_first, nPrevClusters, emptyIndicator] = convertToSurfaceDescriptor(X, lenCombs, layerID, nClusters, n2Clusters, fileForVisualizationPrevLayer{layerID-1},  ...
                       offsetsConventional, depthStep, fieldCenter, cluster1Centres, downsamplingScheme, clusterCurDepths);  % this is a first layer descriptor
                   
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
    
    [coverage] = recomputeCoverage(statistics, outputCoords, outputScales, outputFrames, coverage, X, 1, dx, dy,...
                                    list_depth, lenF, nPrevClusters + 1, lenCombs, is_sparse, triples);
                                
%     coverage(1:kk-1) = ones(1:(kk-1), 1)*(max(coverage)+10);
%     coverageOut(1:kk-1) = ones(1:(kk-1), 1)*(max(coverage)+10);

 
    % sort parts according to their coverage
    [coverage, inds] = sort(coverage, 'descend');
    X = X(inds, :);
    X_first = X_first(inds, :);   
    X = int16(X);
    stat = gpuArray(stat);
    
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
        tic
        

        try
            inds = arrayfun(@isMemberMy, stat(:,1), stat(:,2), stat(:,3));
            elPositions = outputCoords(inds, :);
            scales = outputScales(inds, :);
            frames = outputFrames(inds, :);
            
        catch exception
            
            inds = ismember(stat, curPart, 'rows');
            elPositions = outputCoords(inds, :);
            scales = outputScales(inds, :);
            frames = outputFrames(inds, :);
        end
        
        toc;
      
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
                bw = zeros(r,c,3);

                % the right image is now open
                for k = 1:lenE  % fill all positions for this image
                    
                    centre = double([elPositions(ind(k), 2), elPositions(ind(k), 3)]);
                    frame = reshape(frames(k,:), [3,3]);
                    curScale = double(scales(k, 3));
                    vectX = frame(1:2,2)' * dx;
                    vectY = frame(1:2,3)' * dy;

                    points = round([centre - vectX + vectY; centre + vectX + vectY;  ...
                          centre + vectX - vectY; centre - vectX - vectY]);

                    bw(:,:,3) = poly2mask(points(:,1), points(:,2), r, c);
                    I(bw == 1) = 1;

%                     x = elPositions(ind(k), 2);
%                     y = elPositions(ind(k), 3);
%                     
%                     yMin = max(y-dy, 1);
%                     yMax = min(y+dy, r);
%                     xMin = max(x-dx, 1);
%                     xMax = min(x+dx, c);      
%                     
%                     I(yMin:yMax, xMin:xMax, 3) = 0; 
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
        fullRecFreq = 5;
        
        if mod(i, fullRecFreq) == 0 || (coverage(i+1)<= lenSelected) % full recompute
            lenRec = lenCombs - i;
        elseif mod(i, 2) == 0
            lenRec = 50;
        end

        [coverage] = recomputeCoverage(statistics, outputCoords, outputScales, outputFrames, coverage, X, i+1, dx, dy,...
                                        list_depth, lenF, nPrevClusters + 1, lenRec, is_sparse, triples);                                  
        
        % sort parts according to their coverage (sort only the rest of the array)
        [coverage(i+1:end), inds] = sort(coverage(i+1:end), 'descend');
        inds = inds+i;
        X(i+1:end,:) = X(inds, :);
        X_first(i+1:end,:) = X_first(inds, :);
        
        if (mod(i, fullRecFreq) == 0) && (coverage(i+1)<= lenSelected)  % maximal remaining coverage is less than we expect
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
        
        
    reset(gpuDevice(1));        
end
    





