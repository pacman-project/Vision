% this function returns displacements for recognition of each layer after 4

function [indsXLeft, indsYLeft, indsXRight, indsYRight] = getDisplacements(layerID, centerX, centerY, displ7, displ5, displ3, smallOffset)
    

    if layerID == 5
        % take a square of size 2*smallOffset
        
        len = 3;
        indsXLeft = zeros(1,len^2);
        indsYLeft = zeros(1,len^2);
        indsXRight = zeros(1,len^2);
        indsYRight = zeros(1,len^2);
        
        centerXLeft = centerX - displ5;
        centerYLeft = centerY;
        centerXRight = centerX + displ5;
        centerYRight = centerY;
        
        range = len - 1;
        num = 1;
        
        for i = -range:smallOffset:range      % should be [-4, -2, 0, 2, 4]  % x - direction
            for j = -range:smallOffset:range  % should be [-4, -2, 0, 2, 4]  % y - direction
                
                indsXLeft(num)  = centerXLeft + i;
                indsYLeft(num)  = centerYLeft + j;
                indsXRight(num) = centerXRight + i;
                indsYRight(num) = centerYRight + j;
                num = num + 1;
                
            end
        end
    end

    if layerID == 6
        % take a square of size 2*smallOffset
        
        len = 3;
        indsXLeft  = zeros(1, len^2);
        indsYLeft  = zeros(1, len^2);
        indsXRight = zeros(1, len^2);
        indsYRight = zeros(1, len^2);
        
        centerXLeft = centerX;
        centerYLeft = centerY - displ5;
        centerXRight = centerX;
        centerYRight = centerY + displ5;
        
        range = len - 1;
        num = 1;
        
        for i = -range:smallOffset:range      % should be [-4, -2, 0, 2, 4]  % x - direction
            for j = -range:smallOffset:range  % should be [-4, -2, 0, 2, 4]  % y - direction
                
                indsXLeft(num)  = centerXLeft + i;
                indsYLeft(num)  = centerYLeft + j;
                indsXRight(num) = centerXRight + i;
                indsYRight(num) = centerYRight + j;
                num = num + 1;
                
            end
        end
    end
         
end


