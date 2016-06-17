function [ positions ] = calculatePooledPositions( precisePositions, poolFactor, poolDim, stride )
     tempPositions = floor((double(precisePositions)-1) / stride);
     tempPositions = floor(tempPositions / (poolDim^poolFactor));
     positions = int32(tempPositions + 1);
end