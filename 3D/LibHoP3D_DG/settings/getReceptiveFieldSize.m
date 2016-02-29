% returns radius of the receptive field

function [receptiveFieldRad, offsetsConventional] = getReceptiveFieldSize(dataSetNumber)

    if dataSetNumber < 5   % depth data
        receptiveFieldRad{1} = 2;
        receptiveFieldRad{2} = 7;
        
    elseif dataSetNumber == 5
        
        % this is RADIUS of the receptive field 
        receptiveFieldRad{1} = 0.0033;
        receptiveFieldRad{2} = 0.0033;
        receptiveFieldRad{3} = 0.01;
        receptiveFieldRad{4} = 0.01;
        receptiveFieldRad{5} = 0.03;
        receptiveFieldRad{6} = 0.03;    
        receptiveFieldRad{7} = 0.09;
        receptiveFieldRad{8} = 0.09;
        receptiveFieldRad{9} = 0.27;
        receptiveFieldRad{10} = 0.27;
        
        offsetsConventional{3} = 0.0066;
        offsetsConventional{4} = 0.0066;
        offsetsConventional{5} = 0.02;
        offsetsConventional{6} = 0.02;
        offsetsConventional{7} = 0.06;
        offsetsConventional{8} = 0.06;
        offsetsConventional{9} = 0.18;
        offsetsConventional{10} = 0.18;
    end

end

