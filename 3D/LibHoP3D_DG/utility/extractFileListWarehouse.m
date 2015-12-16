% this is to extract filelist fromt the Washington dataset

function [list_input, categoryID, modelID, lenF] = extractFileListWarehouse(inputPath, is_subset, subsetPercent)
    

    disp('extracting file list...');
    list_input = {};
    categoryID = [];
    modelID = [];

    %  extract a list of subfolders
    subfolders = list_of_subfolders(inputPath);
    % first two elements are garbage
    subfolders(1:2) = [];
    lenS = length(subfolders);

    for i = 1:lenS % for each category
        str = [inputPath, subfolders{i}, '/models/'];

        list_file = load_filelist(str);
        fileD = list_file';
        lenD_cur = length(fileD); 

        if (is_subset == true)  % select images randomly

            if (subsetPercent > 1) || (subsetPercent < 0)
                disp('Error: subset percentage is invalid');
            end

            subset_len = round(lenD_cur * subsetPercent);
            % select the subset
            nums = randperm(lenD_cur);
            nums = nums(1:subset_len);
            list_input = [list_input; fileD(nums)];
            
        else
            list_input = [list_input; fileD];
            categoryID = [categoryID; i*ones(lenD_cur, 1)];
            modelID = [modelID; (1:lenD_cur)'];
        end

    end


    lenF = length(list_input);
end

