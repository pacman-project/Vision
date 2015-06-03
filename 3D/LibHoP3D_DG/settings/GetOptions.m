% options used in learning/inferece procedures

function options = GetOptions(filtOptions)

    options.is2D = false;
    options.methodId = 1;  % paraboloid fitting
    
    options.minValue = filtOptions.minDepth;
    options.maxValue = filtOptions.maxDepth;
    options.step = 2;
    options.planarPatch = 0.04;  % average angular error in degrees
    
    options.maxAngleToZaxis = 65;  % in degrees
    
    % numbe of points to robustly establish the surface
    options.ignoreThresh{1} = 4;
    options.ignoreThresh{2} = 10;
    options.ignoreThresh{3} = 13;
    options.ignoreThresh{4} = 14;
    
    options.weight = 0.5;
    
end

