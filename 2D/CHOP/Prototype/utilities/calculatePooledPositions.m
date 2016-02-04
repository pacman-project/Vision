function [ positions ] = calculatePooledPositions( precisePositions, levelItr, options )
     if strcmp(options.filterType, 'gabor')
          stride = options.gabor.stride;
          inhibitionRadius = options.gabor.inhibitionRadius;
     else
          stride = options.auto.stride;
          inhibitionRadius = options.auto.inhibitionRadius;
     end
     tempPositions = floor((precisePositions-1) / (stride * (inhibitionRadius+1)));
     tempPositions = floor(tempPositions / ((options.poolDim)^(levelItr-2)));
     positions = tempPositions + 1;
end