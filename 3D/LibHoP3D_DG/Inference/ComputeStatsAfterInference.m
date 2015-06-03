% This function computes number of detected elements of each layer after inference procedure 

function [ tableEls ] = ComputeStatsAfterInference( list_El, nPrevClusters)

    lenF = length(list_El);
    tableEls = zeros(1, nPrevClusters);
    
    for i = 1:lenF
        marks = imread(list_El{i});
        
        marks = marks(marks>0);
        
        len = length(marks);
        
        for j = 1:len
            tableEls(marks(j)) = tableEls(marks(j)) + 1;
        end
        
    end


end

