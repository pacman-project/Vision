function [meargeThresh] = defineMeargeThreshes(distanceIndex, dataSetNumber)

if distanceIndex == 1 
    
    if dataSetNumber == 2 || dataSetNumber == 4
        meargeThresh{3} = 0.018;  
        meargeThresh{4} = 0.02; 
        meargeThresh{5} = 0.03; 
        meargeThresh{6} = 0.04; 
        meargeThresh{7} = 0.05; 
        meargeThresh{8} = 0.06; 
        
    elseif dataSetNumber == 1 || dataSetNumber == 3
        
        meargeThresh{3} = 0.1;  
        meargeThresh{4} = 0.2; 
        meargeThresh{5} = 5.0;
        meargeThresh{6} = 9.0; 
        meargeThresh{7} = 25.0; 
        meargeThresh{8} = 25.0;
        
    elseif dataSetNumber == 5
        
        meargeThresh{3} = 0.1;  
        meargeThresh{4} = 0.2; 
        meargeThresh{5} = 0.3;
        meargeThresh{6} = 0.4; 
        meargeThresh{7} = 0.5; 
        meargeThresh{8} = 0.6;
        
    end
    
    
    
end

