% this function returns displacements for recognition of each layer after 4

% elementType may be a circle or a square

function [indsXLeft, indsYLeft, indsXRight, indsYRight] = getDisplacements(layerID, centerX, centerY, displ3, indsX, indsY)
 
    % here we shift indexes to the right central position
    
    if mod(layerID,2) == 1  % 3,5, etc.
        
        indsXLeft = indsX + centerX - displ3; 
        indsYLeft = indsY + centerY;
        indsXRight = indsX + centerX + displ3;
        indsYRight = indsY + centerY;
        
    else
        
        indsXLeft = indsX + centerX; 
        indsYLeft = indsY + centerY - displ3;
        indsXRight = indsX + centerX;
        indsYRight = indsY + centerY + displ3;
        
    end

         
end


