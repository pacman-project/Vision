function [I, Ix, Iy, mask, r, c, is_successfull] = preliminaryProcessing(I, isErrosion, discSize, isY, isTrim, dxKernel, sigmaKernelSize, sigma)

    I = double(I);

    % create a mask
    [r,c] = size(I);
    mask = zeros(r,c);
    mask(I > 0) = 1;

    % check if the mask is empty or not
    maxM = max(max(mask));
    if maxM == 0 
        is_successfull = false;
    else
        if isTrim
            [I, mask] = trimImage(I, mask); % trim the image
            [r,c] = size(I);
        end
        if isErrosion
            se = strel('disk', discSize);  % errosion of the mask
            mask = imerode(mask, se);  
        end

        % compute filter responces
        hx = dxKernel;       
        h = fspecial('gaussian', sigmaKernelSize, sigma);
        IG = imfilter(I, h); % gaussian
        % now computing Gaussian derivatives
        Ix = imfilter(IG, hx);
        Ix = Ix.*mask; % to avoid high on the boundary
        if isY
            hy = hx';
            Iy = imfilter(IG, hy); % derivative in y direction
            Iy = Iy.*mask;
        else
            Iy = [];
        end
        I = IG;
        is_successfull = true;
    end
end