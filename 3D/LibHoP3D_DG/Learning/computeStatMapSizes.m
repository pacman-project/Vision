function [sizeXY, sizeZ, sizeAngle, centreXY, centreZ, centreAngle] = computeStatMapSizes(statMapProperties, offsetConventional)

    xyzStep = statMapProperties.xyzStep;
    vectStep = statMapProperties.vectStep;
    multMax = max(statMapProperties.multX, statMapProperties.multY);
    offsetMaximalAngle = statMapProperties.offsetMaximalAngle;  
    
    sizeXY   = ceil(2*offsetConventional*multMax/xyzStep);
    sizeZ = ceil(2*offsetConventional/xyzStep);
    sizeAngle = round(2*offsetMaximalAngle/vectStep);
    
    if mod(sizeXY, 2) == 0
        sizeXY = sizeXY+1;
    end
    if mod(sizeZ, 2) == 0
        sizeZ = sizeZ+1;
    end
    if mod(sizeAngle, 2) == 0
        sizeAngle = sizeAngle + 1;
    end
    centreXY   = ceil(sizeXY/2);
    centreZ   =  ceil(sizeZ/2);
    centreAngle = ceil(sizeAngle/2);
end

