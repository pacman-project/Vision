% this is the script for greedy part selection for the layer 3.

function [partsOut, coverageOut, lenOut] = partSelectionNew(nClusters, n2Clusters, statistics3LayerSieved, statistics3LayerAggregated, ...
                             dataSetNumber, fieldSize, list_depth, lenF, meargeThresh, ...
                             iterations, layerID, fileForVisualizationPrevLayer)
                         
    %   output variables
    coverageOut = zeros(iterations, 1);
    partsOut = zeros(iterations,3);
    lenOut = 0;
                           
    halfFieldSize = floor(fieldSize/2);    %         for example fieldSize = [17, 5, 71];
    dx = halfFieldSize(1);
    dy = halfFieldSize(2);
    
    load(statistics3LayerSieved);          %       'statistics', 'cluster3Depths', 'outputCoords'
    load(statistics3LayerAggregated);      %       'X' ,'frequencies', 'curTS', 'triples'
    lenStat =  size(statistics, 1);
    lenCombs = size(X, 1);
   
    load(fileForVisualizationPrevLayer);
    downsamplingScheme = 3;

    % apply conversion of parts to actuall surface descriptors
    
    if layerID < 9
        triple8OutDepth = [];
    elseif layerID < 8
        triple7OutDepth = [];
    elseif layerID < 7
        triple6OutDepth = [];
    elseif layerID < 6
        triple5OutDepth = [];
    elseif layerID < 5
        triple4OutDepth = [];
    elseif layerID < 4
        triple3OutDepth = [];
    end
    
    % now define nPrevClusters
    if     layerID == 8
        nPrevClusters = size(triple7OutDepth, 1);
    elseif layerID == 7
        nPrevClusters = size(triple6OutDepth, 1);
    elseif layerID == 6
        nPrevClusters = size(triple5OutDepth, 1);
    elseif layerID == 5
        nPrevClusters = size(triple4OutDepth, 1);
    elseif layerID == 4
        nPrevClusters = size(triple3OutDepth, 1);
    end
     

    [positions, elements] = partMeanReconstruction(layerID, partID, fieldCenter, triple8OutDepth, triple7OutDepth, triple6OutDepth, triple5OutDepth, triple4OutDepth, ...
                                                            triple3OutDepth, displ3, displ5, displ7, nClusters);
                                                        
    [X_first, d] = partToSurfaceDescriptor(elements, positions, nClusters, cluster1Centres, downsamplingScheme);
    
 

    coverage = zeros(lenCombs, 1);   % just something to initialize
    coverage = double(coverage);
    
    [coverage] = recomputeCoverage(statistics, outputCoords, coverage, X, 1, dx, dy,...
                                    list_depth, lenF, nPrevClusters, lenCombs, layerID);


    
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
        
        % merge the selected elements with all elements with distance
        % less than meargeThresh

        distance = Integral_distances(X_first, X_first(i,:), lenCombs, 1, false, false);          
        similarInds = find(distance < meargeThresh);
        curPart = [];
        for j = 1:length(similarInds)
            % extract the element X(similarInds(j), :)
            curPart = [curPart; X(similarInds(j), :)]; 
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
            
            % extract all related to the image
            % 1) find the current position in the elPositions
            ind = find(firstColumn == j);
            lenE = length(ind);
            
            if lenE > 0
                % open the image
                I = imread(list_depth{j});

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
            
        end
        
        coverageOut(i) = coverage(i);
        partsOut(i, :) = X(i,:);
        lenOut = lenOut + 1;
        
 
        %------------------------------------------------------------------
        % update coverage of the following parts
        lenRec = 25;
        
        if mod(i, 3) == 0
            lenRec = 50;
        end
        if mod(i, 10) == 0  % full recompute
            lenRec = lenCombs - i;
        end
        [coverage] = recomputeCoverage(statistics, outputCoords, coverage, X, i+1, dx, dy,...
                                        list_depth, lenF, nPrevClusters, lenRec, layerID);

        % sort parts according to their coverage (sort only the rest of the array)
        [coverage(i+1:end), inds] = sort(coverage(i+1:end), 'descend');
        inds = inds+i;
        X(i+1:end,:) = X(inds, :);
        X_first(i+1:end,:) = X_first(inds, :);
            

    end
        
        
end
    





