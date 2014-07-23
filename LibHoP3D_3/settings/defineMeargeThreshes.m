function [ meargeThresh3, meargeThresh4,meargeThresh5, meargeThresh6, meargeThresh7, meargeThresh8 ] = defineMeargeThreshes(distanceIndex)

if distanceIndex == 1
    meargeThresh3 = sqrt(2) - 0.01;  % sqrt(2) are not mearged
    meargeThresh4 = sqrt(4) - 0.01;
    meargeThresh5 = sqrt(15) - 0.01;
    meargeThresh6 = sqrt(33) - 0.01;
    meargeThresh7 = sqrt(68) - 0.01;
    meargeThresh8 = sqrt(120) - 0.01;
end

