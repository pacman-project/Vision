
function receptiveField = getReceptiveFieldSize(dataSetNumber)

    if dataSetNumber < 5   % depth data
        receptiveField{1} = 2;
        receptiveField{2} = 7;
        
    elseif dataSetNumber == 5  

        receptiveField{1} = 0.15;
        receptiveField{2} = 0.30;
    end

end

