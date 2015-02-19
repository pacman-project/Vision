% the function adds zero boundaries to the image

function [Iout] = addZeroBoundaries(I, boundSize)

    [rm, cm] = size(I);
    if strcmp(class(I), 'double')
        I1 = double(zeros(rm + 2 * boundSize, cm + 2 * boundSize));
    elseif strcmp(class(I), 'uint16')
        I1 = uint16(zeros(rm + 2 * boundSize, cm + 2 * boundSize));
    end

    I1(boundSize+1:end-boundSize, boundSize+1:end - boundSize) = I;
    Iout = I1; 

end

