% if mask is to be created call with mask = []

function [I, Ix, Iy, mask, r, c, is_successfull] = preliminaryProcessing(I, mask, isErrosion, discSize, isX, isY, ...
                    isTrim, dxKernel, sigmaKernelSize, sigma, is_guided, r_guided, eps, is_mask_extended, maxExtThresh1, maxExtThresh2)

    I = double(I);
    [r,c] = size(I);
    is_successfull = true;
        
    % create a mask
    
    if isempty(mask)
        mask = zeros(r,c);
        mask(I > 0) = 1;
    end
    
    % extend mask to regions with no data (for Washington DS only)
    if is_mask_extended
        mask1 = ones(r,c);
        mask1(I < maxExtThresh1 & mask == 1) = 0;  % no data at these points
        mask1(I > maxExtThresh2 & mask == 1) = 0;  % no data at these points
        mask(mask1 == 0) = 0;
    end

    % check if the mask is empty or not
    maxM = max(max(mask));
    if maxM == 0 
        is_successfull = false;
    else
        if isTrim
            [I, mask] = trimImage(I, mask); % trim the image
            [r,c] = size(I);
            
            if min(r,c) < 30
                is_successfull = false;
                Ix = [];
                Iy = [];
            end
        end
        
        if is_successfull
            if isErrosion
                se = strel('disk', discSize);  % errosion of the mask
                mask = imerode(mask, se);  
            end
            if is_guided  % apply guided image filter
                IG = guidedfilter(mask, I, r_guided, eps);
            else % apply gaussian filter
                h = fspecial('gaussian', sigmaKernelSize, sigma);
                IG = imfilter(I, h, 'replicate', 'same'); % gaussian
            end

            if isX
                hx = dxKernel; 
                Ix = imfilter(IG, hx, 'replicate', 'same');
                Ix = Ix.*mask; % to avoid high on the boundary
            else
                Ix = [];
            end

            if isY
                hy = hx';
                Iy = imfilter(IG, hy, 'replicate', 'same'); % derivative in y direction
                Iy = Iy.*mask;
            else
                Iy = [];
            end


            I = IG;
        end
    end
end