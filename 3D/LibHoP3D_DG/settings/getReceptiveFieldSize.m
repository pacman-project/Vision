% returns radius of the receptive field

function receptiveField = getReceptiveFieldSize(dataSetNumber)

    if dataSetNumber < 5   % depth data
        receptiveField{1} = 2;
        receptiveField{2} = 7;
        
    elseif dataSetNumber == 5  

        receptiveField{1} = 0.005;
        receptiveField{2} = 0.005;
        receptiveField{3} = 0.015;
        receptiveField{4} = 0.015;
        receptiveField{5} = 0.045;
        receptiveField{6} = 0.045;
    end

end

