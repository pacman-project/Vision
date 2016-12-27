function [ rfSize, pixelSize ] = getRFSize( options, levelItr )
    if strcmp(options.filterType, 'gabor')
        stride = options.gabor.stride;
    else
        stride = options.auto.stride;
    end
   poolFactor = numel(setdiff(2:(levelItr-1), options.noPoolingLayers));
   pixelSize = stride * round(options.poolDim ^ poolFactor);
   rfSize =  options.receptiveFieldSizes(levelItr-1)  * pixelSize;
   rfSize = [rfSize, rfSize];
end