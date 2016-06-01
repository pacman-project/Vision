function [ positions ] = calculatePooledPositions( precisePositions, poolFactor, poolDim, stride )
     tempPositions = floor((precisePositions-1) / stride);
     tempPositions = floor(tempPositions / (poolDim^poolFactor));
     positions = tempPositions + 1;
end