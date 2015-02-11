% if mask is to be created call with mask = []
% isX_FB - do we need to compute forward and backward difference

% filterType may be:
% 0 - gaussian smoothing
% 1 - guided
% 2 - no filter

function [I, Ix, Iy, mask, r, c, is_successfull] = preliminaryProcessing(I, mask, isErrosion, discSize, isX, isY, isX_FB, ...
                    isTrim, dxKernel, sigmaKernelSize, sigma, filterType, r_guided, eps, is_mask_extended, maxExtThresh1, maxExtThresh2, ...
                    dyKernelTop, dyKernelBottom, dxKernelBack, dxKernelForward)

    I = double(I);
    [r,c] = size(I);
    is_successfull = true;
        
    % create a mask
    
    if isempty(mask)
        mask = zeros(r,c);
        mask(I > 0) = 1;
    end
    
    % extend mask to regions with no data
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
            
            if filterType == 1  % apply guided image filter
                IG = guidedfilter(mask, I, r_guided, eps);
            elseif filterType == 0 % apply gaussian filter
                h = fspecial('gaussian', sigmaKernelSize, sigma);
                IG = imfilter(I, h, 'replicate', 'same'); % gaussian
            end
            
            if filterType ~= 2 % if there is filtering
                I = IG;
            end

            if isX
                hx = dxKernel; 
                Ix = imfilter(IG, hx, 'replicate', 'same');
                Ix = Ix.*mask; % to avoid high on the boundary
                
                if isX_FB % compute forward and backward differences
                    hx = dxKernelBack; 
                    IxB = imfilter(IG, hx, 'replicate', 'same');
                    IxB = IxB.*mask; % to avoid high on the boundary
                    
                    hx = dxKernelForward; 
                    IxF = imfilter(IG, hx, 'replicate', 'same');
                    IxF = IxF.*mask; % to avoid high on the boundary
                    
                    Ix(:,:,2) = IxB;
                    Ix(:,:,3) = IxF;
                    
                end
            else
                Ix = [];
            end

            if isY
                % middle, forward and backward differences
                
                hy = dxKernel';
                Iy = imfilter(IG, hy, 'replicate', 'same'); % derivative in y direction
                Iy = Iy.*mask;
                
                hy = dyKernelTop;
                IyT = imfilter(IG, hy, 'replicate', 'same'); % derivative in y direction
                IyT = IyT.*mask;
                
                hy = dyKernelBottom;
                IyB = imfilter(IG, hy, 'replicate', 'same'); % derivative in y direction
                IyB = IyB.*mask;
                
                Iy(:,:,2) = IyT;
                Iy(:,:,3) = IyB;
            else
                Iy = [];
            end


        end
    end
end