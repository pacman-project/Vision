% this is the script for greedy part selection for the layer 3.

function [partsOut, coverageOut, lenOut] = partSelectionNew(nClusters, n2Clusters, statistics3LayerSieved, statistics3LayerAggregated, ...
                             dataSetNumber, fieldSize, list_depth, lenF, abstractionLevel, abstractionTable3Layer)
                           
    halfFieldSize = floor(fieldSize/2);    %         for example fieldSize = [13, 5, 71];
    dx = halfFieldSize(1);
    dy = halfFieldSize(2);
    
    load(statistics3LayerSieved);          %       'statistics', 'cluster3Depths', 'outputCoords'
    load(statistics3LayerAggregated);      %       'X' ,'frequencies', 'curTS', 'triples'
    
    lenStat =  size(statistics, 1);
    lenCombs = size(X, 1);
    
    if abstractionLevel == 1 || abstractionLevel == 2 
        is_abstraction = true;
    else
        is_abstraction = false;
    end
    
    if is_abstraction
        [abstractionTable] = computeAbstractionTable(X, nClusters, n2Clusters, abstractionLevel);
        save(abstractionTable3Layer, 'abstractionTable');  % save the table for the unsorted X
    else
        abstractionTable = [];
    end
 
%   output variables
    coverageOut = zeros(1,lenCombs);
    partsOut = zeros(lenCombs,3);
    lenOut = 0;
    
    recompute = ones(1, n2Clusters);
    
    coverage = zeros(lenCombs, 1);   % just something to initialize
    coverage = double(coverage);
    
    [coverage] = recomputeCoverage(statistics, abstractionTable, outputCoords, coverage, X, 1, recompute, dx, dy, list_depth, lenF, n2Clusters);
    recompute = recompute * 0;
    
    % sort parts according to their coverage
    [coverage, inds] = sort(coverage, 'descend');
    X = X(inds, :);
    
    if is_abstraction   
        % sort the abstraction table
        for j = 1:lenCombs
            % replace pointers to the elements in the abstractionTable
            if abstractionTable(j) ~= 0 % there is a pointer
                newPos = find(inds == abstractionTable(j));
                abstractionTable(j) = newPos(1);
            end
        end
        abstractionTable = abstractionTable(inds);
    end
    
    for i = 1:800   % for each part 

        str = ['element - ', num2str(i)];
        disp(str);
             
        elPositions = []; % positions of this part in all images
        curElPosition = 1;
         
        curPart = X(i,:); % take the part with the largest coverage
        
        if is_abstraction % find the positions of the similar elements in images
            similarInds = find(abstractionTable == i); % all parts that refer to X(i,:)
            
            for j = 1:length(similarInds)
                % extract the element X(similarInds(j), :)
                curPart = [curPart; X(similarInds(j), :)];    
            end
            
        end
        
        for j = 1:size(curPart,1) % update recompute list
            for k = 1:size(curPart,2)
                recompute(curPart(j,k)) = 1;
            end
        end
        
        numEls = size(curPart, 1);
        % find the parts's positions in all images
        parfor j = 1:lenStat
            
            stat = statistics(j,:);
            line = [stat(2), stat(1), stat(4)];
            % match curPart and line
            for k = 1:numEls
                if isequal(line, curPart(k,:))
                    elPositions = [elPositions ; outputCoords(j, :)];
                    break;
                end
            end
        end
        
        % project this part to the images
        % curElPosition is already sorted by image number
        firstColumn = elPositions(:,1);
        
        parfor j = 1:lenF
            % open the image
            I = imread(list_depth{j});
            
            % extract all reated to the image
            % 1) find the current position in the elPositions
            ind = find(firstColumn == j);
            lenE = length(ind);

            % the right image is now open
            for k = 1:lenE  % fill all positions for this image
            
                x = elPositions(ind(k), 2);
                y = elPositions(ind(k), 3);
                I(y-dy:y+dy ,x-dx:x+dx, 3) = I(y-dy:y+dy ,x-dx:x+dx, 3) * 0;
                % I(y-dy:y+dy ,x-dx:x+dx, 2) = I(y-dy:y+dy ,x-dx:x+dx, 2) * 0;
            end

            % save the last image
            I = uint16(I);
            imwrite(I, list_depth{j}, 'png');
            
%             if mod(j,10) == 0
%                 disp('image');
%                 j
%             end
        end
        
        coverageOut(i) = coverage(i);
        partsOut(i, :) = X(i,:);
        lenOut = lenOut + 1;
        
        % when part is selected: destroy links to this group in the
        % abstraction table
        groupA = find(abstractionTable == i);  % define a group of elements that refer to X(i,:)
        % destroy all links to the group elements
        lenGr = length(groupA);
        for jj = 1:lenGr
            subGroup = abstractionTable == groupA(jj);  % all parts which refer to groupA(jj)
            abstractionTable(subGroup) = 0;
        end
        
        %------------------------------------------------------------------
        % update coverage of the following parts
        
        nextPart = X(i+1,:);
        
        if recompute(nextPart(1)) == 1 || recompute(nextPart(2)) == 1 || recompute(nextPart(3)) == 1 % updates required
            
            recompute = ones(1, n2Clusters);
            [coverage] = recomputeCoverage(statistics, abstractionTable, outputCoords, coverage, X, i+1, recompute, dx, dy, list_depth, lenF, n2Clusters);
            recompute = recompute * 0;
            
            % sort parts according to their coverage-----------------------
            [coverage, inds] = sort(coverage, 'descend');
            X = X(inds, :);
            
            if is_abstraction  
                % sort the abstraction table
                for j = 1:lenCombs
                    % replace pointers to the elements in the abstractionTable
                    if abstractionTable(j) ~= 0 % there is a pointer
                        newPos = find(inds == abstractionTable(j));
                        abstractionTable(j) = newPos(1);
                    end
                end
                abstractionTable = abstractionTable(inds);
            end
            %--------------------------------------------------------------
            
            a = 2;

        end
        
        
    end
    

end



