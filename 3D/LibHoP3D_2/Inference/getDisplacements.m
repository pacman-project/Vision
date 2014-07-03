% this function returns displacements for recognition of each layer after 4

function [indsXLeft, indsYLeft, indsXRight, indsYRight] = getDisplacements(layerID, centerX, centerY, displ7, displ5, displ3, hDispl)
    

    if layerID == 5
         
        indsXLeft  = [centerX - displ5, centerX - displ5 - hDispl, centerX - displ5 + hDispl, centerX - displ5, centerX - displ5 - hDispl, centerX - displ5 + hDispl, centerX - displ5, centerX - displ5 - hDispl, centerX - displ5 + hDispl];
        indsYLeft  = [centerY, centerY, centerY, centerY - hDispl, centerY - hDispl, centerY - hDispl, centerY + hDispl, centerY + hDispl, centerY + hDispl];
        
        indsXRight = [centerX + displ5, centerX + displ5 + hDispl, centerX + displ5 - hDispl, centerX + displ5, centerX + displ5 + hDispl, centerX + displ5 - hDispl, centerX + displ5, centerX + displ5 + hDispl, centerX + displ5 - hDispl];
        indsYRight = [centerY, centerY, centerY, centerY - hDispl, centerY - hDispl, centerY - hDispl, centerY + hDispl, centerY + hDispl, centerY + hDispl];

    end

    if layerID == 6   % TODO
         
    end

end

