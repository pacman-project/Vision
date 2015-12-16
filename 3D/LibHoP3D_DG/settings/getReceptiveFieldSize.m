% returns radius of the receptive field

function [receptiveFieldRad, offsetsConventional] = getReceptiveFieldSize(dataSetNumber)

    if dataSetNumber < 5   % depth data
        receptiveFieldRad{1} = 2;
        receptiveFieldRad{2} = 7;
        
    elseif dataSetNumber == 5
        
        % this is RADIUS of the receptive field 
        receptiveFieldRad{1} = 0.01;
        receptiveFieldRad{2} = 0.01;
        receptiveFieldRad{3} = 0.03;
        receptiveFieldRad{4} = 0.03;
        receptiveFieldRad{5} = 0.09;
        receptiveFieldRad{6} = 0.09;
        
        offsetsConventional{3} = 0.02;
        offsetsConventional{4} = 0.02;
        offsetsConventional{5} = 0.06;
        offsetsConventional{6} = 0.06;
        offsetsConventional{5} = 0.2;
        offsetsConventional{6} = 0.2;
    end

end

