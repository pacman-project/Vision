function [quantilesFirst] = defineFirstLayerQuantiles(nClusters, dataSetNumber, is_guided)

    if ~is_guided  % parametres for gaussian filtering
        if nClusters == 9
            if dataSetNumber == 1
                quantilesFirst = [0.01 0.054 0.11 0.22 0.33];
            elseif dataSetNumber == 2
                quantilesFirst = [0.01 0.045 0.075 0.14 0.30];
            elseif dataSetNumber == 3
                quantilesFirst = [0.01 0.050 0.11 0.22 0.37];
            end

        elseif nClusters == 7
            if dataSetNumber == 1
                quantilesFirst = [0.01 0.054 0.13 0.26];  
            elseif dataSetNumber == 2
                quantilesFirst = [0.02 0.06 0.16 0.40];
            elseif dataSetNumber == 3
                quantilesFirst = [0.01 0.060 0.13 0.26];
            end
        end
        
        
    elseif is_guided % parametres for guided image filtering
        if nClusters == 9
            if dataSetNumber == 1
                quantilesFirst = [0.01 0.064 0.11 0.22 0.33];
            elseif dataSetNumber == 2
                quantilesFirst = [0.01 0.045 0.075 0.14 0.30];
            elseif dataSetNumber == 3
                quantilesFirst = [0.01 0.050 0.11 0.22 0.37];
            end

        elseif nClusters == 7
            if dataSetNumber == 1
                quantilesFirst = [0.01 0.064 0.13 0.26];  
            elseif dataSetNumber == 2
                quantilesFirst = [0.02 0.06 0.16 0.40];
            elseif dataSetNumber == 3
                quantilesFirst = [0.01 0.060 0.13 0.26];
            end
        end
    end
    
    

end

