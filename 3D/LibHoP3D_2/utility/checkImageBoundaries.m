% checks if there are some indexes outside image borders

function [indsXLeft, indsYLeft, indsXRight, indsYRight] = checkImageBoundaries(indsXLeft, indsYLeft, indsXRight, indsYRight, r, c)

    lenL = length(indsXLeft);
    indsTrueL = true(1, lenL);
    
    % check indexes
    indsTrueL(indsXLeft <= 0) = false;
    indsTrueL(indsYLeft <= 0) = false;
    
    % filter out
    indsXLeft = indsXLeft(indsTrueL);
    indsYLeft = indsYLeft(indsTrueL);
    
    lenR = length(indsXRight);
    indsTrueR = true(1, lenR);
    
    % check indexes
    indsTrueR(indsXRight > c) = false;
    indsTrueR(indsYRight > r) = false;
    
    % filter out
    indsXRight = indsXRight(indsTrueR);
    indsYRight = indsYRight(indsTrueR);
    
end

