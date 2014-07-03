% this is to extract filelist fromt the Washington dataset

function [list_el, model_id, category_id, lenF] = extractFileListWashingtonForClassification(strRootE, is_subset, subsetPercent)

    disp('extracting file list...');
    list_el = {};
    model_id = [];
    category_id = [];

    %  extract a list of subfolders
    subfolders = list_of_subfolders(strRootE);
    % first two elements are garbage
    subfolders(1:2) = [];
    lenS = length(subfolders);

    parfor i = 1:lenS % for each category
        str = [strRootE, '/', subfolders{i}];
        % select a list of subfolders again
        subSubFolders = list_of_subfolders(str);
        subSubFolders(1:2) = [];
        lenSS = length(subSubFolders);

        for j = 1:lenSS  % for each object

            strS = [str, '/', subSubFolders{j}];
            list_file = load_filelist(strS);
            list_file = list_file';
            lenSSF = length(list_file);

            % select images randomly
            if (is_subset == true)
                if (subsetPercent > 1) || (subsetPercent < 0)
                    disp('Error: subset percentage is invalid');
                    %return;
                end
                subset_len = round(lenSSF * subsetPercent);
                % select the subset
                nums = randperm(lenSSF);
                nums = nums(1:subset_len);
                
                curMID = ones(subset_len, 1) * j;
                curCID = ones(subset_len, 1) * i;
                
                model_id = [model_id; curMID];
                category_id = [category_id; curCID];
                list_el = [list_el; list_file(nums)];

            else
                curMID = ones(lenSSF, 1) * j;
                curCID = ones(lenSSF, 1) * i;
                
                category_id = [category_id; curCID];
                model_id = [model_id; curMID];
                list_el = [list_el; list_file];
            end

        end
        i          
    end
        
    lenF = length(list_el);
    
end

