% this is the file to project all statistics to the images

function [coverageOverall, areaOverall] = ProjectAllStatistics(statisticsSieved, statisticsAggregated, fieldSize, ...
                                        list_depth, list_mask, lenF, dataSetNumber)

    % statistics files   
    
    checkImages(list_depth, lenF);
    
    scaleFactor = uint16(floor(65535 / 81)); 

    load(statisticsSieved);  % should contain variables 'outputCoords' and 'statistics'        
    load(statisticsAggregated);  % 'X' ,'frequencies', 'curTS', 'triples'
    
    halfFieldSize = floor(fieldSize/2);   % for example fieldSize = [17, 5, 71];
    dx = halfFieldSize(1);
    dy = halfFieldSize(2);
    
    coverageOverall = 0;
    coverageOverall = double(coverageOverall);
    areaOverall = 0;
    areaOverall = double(areaOverall);
    
    firstColumn = outputCoords(:, 1);
    
%     counter = 0;
%     imageIds = [];
    
    for j = 1:lenF
        % open the image
        I = imread(list_depth{j});

        % extract all elements reated to the image
        % 1) find the current position in the elPositions
        ind = find(firstColumn == j);
        lenE = length(ind);
        
        I1 = I(:,:,1);
        I2 = I(:,:,2) * 0;

        % the right image is now open
        for k = 1:lenE  % fill all positions for this image
            
            % this fragment can be used to find a particular element in all images
%             centre = statistics(ind(k),1);
%             left = statistics(ind(k),2);
%             right = statistics(ind(k),4);
%             
%             if left == 5 && centre == 41 && right == 68

                x = outputCoords(ind(k), 2);
                y = outputCoords(ind(k), 3);

                a = scaleFactor * uint16(ones(fieldSize(2), fieldSize(1)));

%                 I(y-dy:y+dy ,x-dx:x+dx, 1) = uint16(a * uint16(statistics(ind(k),1)));
%                 I(y-dy:y+dy ,x-dx:x+dx, 2) = uint16(a * uint16(statistics(ind(k),2)));
                I(y-dy:y+dy ,x-dx:x+dx, 3) = 0; %uint16(a * uint16(statistics(ind(k),4)));
                I2(y-dy:y+dy ,x-dx:x+dx)  = ones(fieldSize(2), fieldSize(1)) * 999;
                
%                 counter = counter + 1;
%                 imageIds = [imageIds, j];
                
        end

        % save the last image
        I = uint16(I);
        imwrite(I, list_depth{j}, 'png');
        I2(I1 == 0) = 0;   % area covered by empty cells does not count
        
        cover = length(I2(I2 == 999));
        if dataSetNumber == 1 || dataSetNumber == 3
            area = length(I1(I1 > 0));
        elseif dataSetNumber == 2
            mask = imread(list_mask{j});
            area = length(mask > 0); % area of the object
        end
 
        coverageOverall = coverageOverall + cover;
        areaOverall = areaOverall + area;

        
        if mod(j,10) == 0
            j
        end
        
    end
    
%     counter

end