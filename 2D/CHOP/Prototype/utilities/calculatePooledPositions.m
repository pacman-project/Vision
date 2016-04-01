function [ positions ] = calculatePooledPositions( precisePositions, levelItr, options )
     if strcmp(options.filterType, 'gabor')
          stride = options.gabor.stride;
     else
          stride = options.auto.stride;
     end
     poolFactor = numel(setdiff(2:(levelItr-1), options.noPoolingLayers));
     tempPositions = floor((precisePositions-1) / stride);
     tempPositions = floor(tempPositions / ((options.poolDim)^poolFactor));
     positions = tempPositions + 1;
end