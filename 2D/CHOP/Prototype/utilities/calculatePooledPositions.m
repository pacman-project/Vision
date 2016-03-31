function [ positions ] = calculatePooledPositions( precisePositions, levelItr, options )
     if strcmp(options.filterType, 'gabor')
          stride = options.gabor.stride;
     else
          stride = options.auto.stride;
     end
     tempPositions = floor((precisePositions-1) / stride);
     tempPositions = floor(tempPositions / ((options.poolDim)^(levelItr-2)));
     positions = tempPositions + 1;
end