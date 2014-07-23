function [is_ok] = VisualizeOverlappedPlements(partIds, partLayers, partCentrePositions, partColors, fieldSize, fileForVisualization4Layer)

    if (nargin < 1)  % default parameters
        
        fileForVisualization4Layer = './statistics/fileForVisualization_4_3_9.mat'; 
        statistics1Layer = './statistics/statistics_1_3_9.mat';
        fieldSize = [27, 27, 91];
        
        
        partsIds = [1,2,3];
        partsLayers = [4,3,3];
        partCentrePositions = [14, 11, 45;   14, 18, 48;  16, 21, 46];
        partColors = [1, 0, 0;   0, 1, 0;   0, 0, 1];
        
    end
    
    displ = 6;
    offsetY = [-displ, 0, displ];
    
    load(fileForVisualization4Layer);   % 'triple3OutDepth', 'triple4OutDepth'; are there
    load(statistics1Layer);   %  'cluster1Centres', 'cluster1Lengths', 'thresh', 'nClusters', 'dataSetNumber'
    
    indsEls4 = [1,5,6];
    depthStep = thresh / 4;

    f = figure;
    nEls = length(partsIds);
    
    table = zeros(nClusters, nClusters);

    for i = 1:nClusters
        for j = 1:nClusters
            table(i,j) = compute2elementIndex(i, j, nClusters);
        end
    end
            
    for j = 1:nEls  % for each 4th layer element
        
        if partsLayers(j) == 4
            cur4triple = triple4OutDepth(partsIds(j), :);
            els4 = cur4triple(indsEls4);
            depths4Adder = [cur4triple(4), 0 cur4triple(9)];  % top centre bottom

            for i = 1:3                            % three triples of the third layet
                curEl = triple3OutDepth(els4(i),:);
                indsEl = [1,2,6,7,8,9];
                indsDepths = [3,4,5,10,11,12];

                curXY = curEl(indsEl);
                depths = curEl(indsDepths);

                elements = zeros(1,3);
                for jj = 1:3
                    elements(jj) = table(curXY(2*jj-1), curXY(2*jj));  % part index in range [0, n2Clusters]
                end
                % define positions

                positionLeft   = [partCentrePositions(j,1) - displ, partCentrePositions(j,2) + offsetY(i), partCentrePositions(j,3) + depths(3) + depths4Adder(i)];
                positionCentre = [partCentrePositions(j,1)        , partCentrePositions(j,2) + offsetY(i), partCentrePositions(j,3)             + depths4Adder(i)];
                positionRight =  [partCentrePositions(j,1) + displ, partCentrePositions(j,2) + offsetY(i), partCentrePositions(j,3) + depths(6) + depths4Adder(i)];

                positions = [positionLeft; positionCentre; positionRight]; % left, middle, right
                surfaceVisualizer(fieldSize, positions, elements, nClusters, cluster1Centres, depthStep, partColors(j,:));
                hold on

            end
            
        elseif partsLayers(j) == 3
            curEl = triple3OutDepth(partsIds(j),:);
            indsEl = [1,2,6,7,8,9];
            indsDepths = [3,4,5,10,11,12];

            curXY = curEl(indsEl);
            depths = curEl(indsDepths);

            elements = zeros(1,3);
            for jj = 1:3
                elements(jj) = table(curXY(2*jj-1), curXY(2*jj));  % part index in range [0, n2Clusters]
            end
            % define positions
            partCenter =    [partCentrePositions(j,1),         partCentrePositions(j,2), partCentrePositions(j,3)];
            positionLeft =  [partCentrePositions(j,1) - displ, partCentrePositions(j,2), partCentrePositions(j,3) + depths(3)];
            positionRight = [partCentrePositions(j,1) + displ, partCentrePositions(j,2), partCentrePositions(j,3) + depths(6)];

            positions = [positionLeft; partCenter; positionRight]; % left, middle, right
            surfaceVisualizer(fieldSize, positions, elements, nClusters, cluster1Centres, depthStep, partColors(j,:));
            hold on
        end
    end
    
    
    is_ok = true;
end

