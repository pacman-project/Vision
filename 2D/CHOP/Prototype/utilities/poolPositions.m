function [ positions ] = poolPositions( positions, poolDim )
     tempPositions = floor(double((positions-1)) / round(poolDim));
     positions = int32(tempPositions + 1);
end