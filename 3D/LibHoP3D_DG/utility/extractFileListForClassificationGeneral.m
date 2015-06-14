% this function extracts a file list (files of parts) for different
% dataSets


function [list_el, model_id, category_id, lenF] = extractFileListForClassificationGeneral(strRootE, is_subset, subsetPercent, dataSetNumber)

    % extracts files 
    if dataSetNumber == 2
        [list_el, model_id, category_id, lenF] = extractFileListWashingtonForClassification(strRootE, is_subset, subsetPercent);   
    elseif dataSetNumber == 4  % Mirela dataSet       
        [list_el, model_id, category_id, lenF] = extractFileListMirelaForClassification(strRootE, is_subset, subsetPercent);  % no masks and no RGBD  
    else
        [list_el, lenF] = extractFileList(strRootE, is_subset, subsetPercent);
        model_id = [];
        category_id = [];
    end

end

