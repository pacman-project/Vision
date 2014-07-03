% this is the script for greedy part selection for the layer 3.

function [partsOut, coverageOut, lenOut] = partSelectionNew(nClusters, n2Clusters, statistics3LayerSieved, statistics3LayerAggregated, ...
                             dataSetNumber, fieldSize, list_depth, lenF, abstractionLevel, abstractionTable3Layer, meargeThresh, ...
                             iterations, layerID, fileForVisualizationPrevLayer)
                         
    %   output variables
    coverageOut = zeros(1,iterations);
    partsOut = zeros(iterations,3);
    lenOut = 0;
                           
    halfFieldSize = floor(fieldSize/2);    %         for example fieldSize = [17, 5, 71];
    dx = halfFieldSize(1);
    dy = halfFieldSize(2);
    
    load(statistics3LayerSieved);          %       'statistics', 'cluster3Depths', 'outputCoords'
    load(statistics3LayerAggregated);      %       'X' ,'frequencies', 'curTS', 'triples'
    lenStat =  size(statistics, 1);
    lenCombs = size(X, 1);
   
    
    if layerID == 3  
        
        table2 = zeros(n2Clusters, 2);
        for i = 1:n2Clusters
            [clusterX, clusterY] = compute2derivatives(i, nClusters);
            table2(i, 1) = clusterX;
            table2(i, 2) = clusterY;
        end
        X_first = Convert3ToFirstLayer(X, lenCombs, table2);  % to express parts in terms of the first layer element
        nPrevClusters = n2Clusters;
        
    elseif layerID == 4
        [X_first, nPrevClusters] = Convert4ToFirstLayer(X, lenCombs, fileForVisualizationPrevLayer); %  'triple3OutDepth'
    elseif layerID == 5
        [X_first, nPrevClusters] = Convert5ToFirstLayer(X, lenCombs, fileForVisualizationPrevLayer); %  'triple3OutDepth', 'triple4OutDepth'
    elseif layerID == 6
        [X_first, nPrevClusters] = Convert6ToFirstLayer(X, lenCombs, fileForVisualizationPrevLayer); %  'triple3OutDepth', 'triple4OutDepth', 'triple5OutDepth'
    end
    
    coverage = zeros(lenCombs, 1);   % just something to initialize
    coverage = double(coverage);
    
    tic;
    [coverage] = recomputeCoverage(statistics, outputCoords, coverage, X, 1, dx, dy,...
                                    list_depth, lenF, nPrevClusters, lenCombs);
    toc

    
    % sort parts according to their coverage
    [coverage, inds] = sort(coverage, 'descend');
    X = X(inds, :);
    X_first = X_first(inds, :);   
    X = int16(X);
    
    for i = 1:iterations   % for each part 
        
        str = ['element - ', num2str(i)];
        disp(str);
        
        if coverage(i) == 0
            return;  % all parts are selected
        end 
        
        if abstractionLevel == 3  % 
            % merge the selected elements with all elements with distance
            % less than meargeThresh
            
            distance = Isodata_distances(X_first, X_first(i,:), lenCombs, 1, false, false);          
            similarInds = find(distance < meargeThresh);
            curPart = [];
            for j = 1:length(similarInds)
                % extract the element X(similarInds(j), :)
                curPart = [curPart; X(similarInds(j), :)]; 
            end 
        end
        
        numEls = size(curPart, 1);
        
        % find the parts's positions in all images 
        elPositions = []; % positions of this part in all images
        indsLine = [2,1,4];
        stat = statistics(:, indsLine);
        

        for k = 1:numEls
            % match this line to all lines in stat
            kk = repmat(curPart(k,:), lenStat, 1);
            a = abs(stat - kk);
            a = sum(a,2);
            inds = a == 0;
            elPositions = [elPositions; outputCoords(inds, :)];
        end

        
        % project this part to the images
        % curElPosition is already sorted by image number
        firstColumn = elPositions(:,1);
        
        parfor j = 1:lenF
            % open the image
            I = imread(list_depth{j});
            
            % extract all related to the image
            % 1) find the current position in the elPositions
            ind = find(firstColumn == j);
            lenE = length(ind);

            % the right image is now open
            for k = 1:lenE  % fill all positions for this image
            
                x = elPositions(ind(k), 2);
                y = elPositions(ind(k), 3);
                I(y-dy:y+dy ,x-dx:x+dx, 3) = zeros(2*dy+1, 2*dx+1); 
            end

            % save the last image
            I = uint16(I);
            imwrite(I, list_depth{j}, 'png');
            
        end
        
        coverageOut(i) = coverage(i);
        partsOut(i, :) = X(i,:);
        lenOut = lenOut + 1;
        
 
        %------------------------------------------------------------------
        % update coverage of the following parts
        lenRec = 10;
        
        if mod(i, 5) == 0
            lenRec = 15;
        end
        if mod(i, 10) == 0
            lenRec = 25;
        end
        if mod(i, 20) == 0  % full recompute
            lenRec = lenCombs - i;
        end
        [coverage] = recomputeCoverage(statistics, outputCoords, coverage, X, i+1, dx, dy,...
                                        list_depth, lenF, nPrevClusters, lenRec);

        % sort parts according to their coverage-----------------------
        [coverage, inds] = sort(coverage, 'descend');
        X = X(inds, :);
        X_first = X_first(inds, :);
            

    end
        
        
end
    





