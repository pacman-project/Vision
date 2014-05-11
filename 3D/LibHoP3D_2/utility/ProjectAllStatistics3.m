% this is the file to project all statistics to the images

function [coverageOverall] = ProjectAllStatistics3(statistics3Sieved, statistics3Aggregated, fieldSize, list_depth, lenF)

    % statistics files   
    
    scaleFactor = uint16(floor(65535 / 81)); 
    
    checkImages(list_depth, lenF);
    
    load(statistics3Sieved);  % should contain variables 'outputCoords' and 'statistics'        
    load(statistics3Aggregated);  % 'X' ,'frequencies', 'curTS', 'triples'
    
    halfFieldSize = floor(fieldSize/2);   % for example fieldSize = [13, 5, 71];
    dx = halfFieldSize(1);
    dy = halfFieldSize(2);
    
    coverageOverall = 0;
    coverageOverall = double(coverageOverall);
    
    firstColumn = outputCoords(:, 1);
    
    for j = 1:lenF
        % open the image
        I = imread(list_depth{j});

        % extract all elements reated to the image
        % 1) find the current position in the elPositions
        ind = find(firstColumn == j);
        lenE = length(ind);
        
        I2 = I(:,:,2) * 0;

        % the right image is now open
        for k = 1:lenE  % fill all positions for this image

            x = outputCoords(ind(k), 2);
            y = outputCoords(ind(k), 3);
            
            a = scaleFactor * uint16(ones(fieldSize(2), fieldSize(1)));
            
            I(y-dy:y+dy ,x-dx:x+dx, 1) = uint16(a * uint16(statistics(ind(k),1)));
            I(y-dy:y+dy ,x-dx:x+dx, 2) = uint16(a * uint16(statistics(ind(k),2)));
            I(y-dy:y+dy ,x-dx:x+dx, 3) = uint16(a * uint16(statistics(ind(k),4)));
            I2(y-dy:y+dy ,x-dx:x+dx)  = ones(fieldSize(2), fieldSize(1)) * 999;
        end

        % save the last image
        I = uint16(I);
        imwrite(I, list_depth{j}, 'png');
        
        cover = length(I2(I2 == 999));
 
        coverageOverall = coverageOverall + cover;
        
        if mod(j,10) == 0
            j
        end
        
    end

end