function [meargeThresh] = defineMeargeThreshes(distanceIndex, dataSetNumber)

if distanceIndex == 1 
    
    if dataSetNumber == 2
        meargeThresh{3} = 0.018;  
        meargeThresh{4} = 0.02; 
        meargeThresh{5} = 0.02; 
        meargeThresh{6} = 0.02; 
        meargeThresh{7} = 0.02; 
        meargeThresh{8} = 0.02; 
        
    elseif dataSetNumber == 1 || dataSetNumber == 3
        
        meargeThresh{3} = 0.1;  
        meargeThresh{4} = 0.2; 
        meargeThresh{5} = 5.0;
        meargeThresh{6} = 9.0; 
        
        
        meargeThresh{7} = 25.0; 
        meargeThresh{8} = 25.0;
    end
end

