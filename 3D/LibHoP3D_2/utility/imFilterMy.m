function [I_filt] = imFilterMy(I, kernel)  % this is for the filter kernel 5x1
    

    kernel = double(kernel);
    [rK, cK] = size(kernel);  % size of the kernel
    if (rK ~= 1) || (cK ~= 5)
        I_filt = imfilter(I, kernel);
    else
        offset = 2;
        I = addZeroBoundaries(I, offset);
        [r,c] = size(I);
        I_filt = zeros(r,c);
        
        
        for i = offset+1:r-offset
            for j = offset+1:c-offset  
                temp = I(i, j-offset:j+offset) .* kernel;
                I_filt(i,j) = sum(temp(:));
            end
        end
        
        I_filt = I_filt(offset+1 : r-offset, offset+1 : c-offset); % remove zero boundaries
    end
end