% the function adds zero boundaries to the image

function [Iout] = addZeroBoundaries(I, boundSize)

    [rm, cm] = size(I);
    I1 = zeros(rm + 2 * boundSize, cm + 2 * boundSize);
    I1(boundSize+1:end-boundSize, boundSize+1:end - boundSize) = I;
    Iout = I1; 

end

