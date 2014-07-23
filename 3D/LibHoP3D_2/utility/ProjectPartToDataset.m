% this is the script for greedy part selection for the layer 3.
% partId should be a triple!


function [listDepthOut, listMaskOut] = ProjectPartToDataset(nClusters, n2Clusters, statistics3LayerSieved, statistics3LayerAggregated, ...
                             fieldSize, list_depth, list_mask, lenF, partID)
                         
    listDepthOut = {};
    listMaskOut = {};
    lenLD = 0;
    halfFieldSize = floor(fieldSize/2);    %         for example fieldSize = [17, 5, 71];
    dx = halfFieldSize(1);
    dy = halfFieldSize(2);
    
    load(statistics3LayerSieved);          %       'statistics', 'cluster3Depths', 'outputCoords'
    load(statistics3LayerAggregated);      %       'X' ,'frequencies', 'curTS', 'triples'
    lenStat =  size(statistics, 1);
    lenCombs = size(X, 1);
   
    curPart = int16(partID);
    % find the parts's positions in all images 
    elPositions = []; % positions of this part in all images
    indsLine = [2,1,4];
    stat = statistics(:, indsLine);
    

    for k = 1:1
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

    for j = 1:lenF

        % 1) find the current position in the elPositions
        ind = find(firstColumn == j);
        lenE = length(ind);

        if lenE > 0
            
            lenLD = lenLD + 1;
            listDepthOut{lenLD} = list_depth{j};
            listMaskOut{lenLD} = list_mask{j};
            % open the image
            I = imread(list_depth{j});

            % the right image is now open
            for k = 1:lenE  % fill all positions for this image

                x = elPositions(ind(k), 2);
                y = elPositions(ind(k), 3);
                I(y-dy:y+dy ,x-dx:x+dx, 3) = zeros(2*dy+1, 2*dx+1); 
            end
            
            % save the last image
            I = I*50;
            
            str = ['D:\Input Data\Trial\',int2str(j), '.png'];
            imwrite(I, str, 'png');
        end

    end
        
end
    





