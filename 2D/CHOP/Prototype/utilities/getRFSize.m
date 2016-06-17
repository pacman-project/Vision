function [ rfSize ] = getRFSize( options, levelItr )
    if strcmp(options.filterType, 'gabor')
        stride = options.gabor.stride;
    else
        stride = options.auto.stride;
    end
   poolFactor = numel(setdiff(2:(levelItr-1), options.noPoolingLayers));
   rfSize = options.receptiveFieldSize * stride * (options.poolDim ^ poolFactor);
   rfSize = [rfSize, rfSize];
end