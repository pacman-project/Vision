function TestDifferentialGeometry()

    I = imread('D:/Input Data/AimShape/1_views_1Scale/1_1_1_1_1_1.png');
    
    I = I(:,:,1);
    
    dxKernel = 0.5 * [-1, 0, 1];
    dyKernel = dxkernel';
    
    % compute first derivative at each point
    
    Ix = imfilter(I, dxKernel);
    Iy = imfilter(I, dyKernel);
    
    % first fundamental form
    
    
    
    
    
    
    
    imtool(I, [0, 2^16]);

end

