% this is to extract filelist fromt the Washington dataset

function [list_el, model_id, category_id, lenF] = extractFileListMirelaForClassification(depthPath, is_subset, subsetPercent)
    
    is_RGB = false;
    is_Mask = false;
    list_el = {};
    
    model_id = [];
    category_id = [];

    disp('extracting file list...');

    if is_RGB
        list_RGB = {};
    else
        list_RGB = [];
    end
    
    if is_Mask
        list_mask = {};
    else
        list_mask = [];
    end

    %  extract a list of subfolders
    subfolders = list_of_subfolders(depthPath);
    % first two elements are garbage
    subfolders(1:2) = [];
    lenS = length(subfolders);

    parfor i = 1:lenS % for each category
        str = [depthPath, '/', subfolders{i}];
        % select a list of subfolders again
        subSubFolders = list_of_subfolders(str);
        subSubFolders(1:2) = [];
        lenSS = length(subSubFolders);

        for j = 1:lenSS  % for each object

            strS = [str, '/', subSubFolders{j}];
            
            if is_RGB
                strSRGB = [strS, '/rgb'];
                list_file = load_filelist(strSRGB);
                fileRGB = list_file';
            end
            
%             strD= [strS, '/depth'];
            list_file = load_filelist(strS);
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

                list_el = [list_el; fileD(nums)];
                
                if is_RGB
                    list_RGB =   [list_RGB;   fileRGB(nums)];
                end
                if is_Mask
%                   list_mask =  [list_mask;  fileMask(nums)];
                end
                
                % define modelID and categoryID
                curMID = ones(subset_len, 1) * j;
                curCID = ones(subset_len, 1) * i;
                
                model_id = [model_id; curMID];
                category_id = [category_id; curCID];
                
            else
                list_el = [list_el; fileD];
                if is_RGB
                    list_RGB =   [list_RGB;   fileRGB];
                end
                if is_Mask
%                 list_mask =  [list_mask;  fileMask];
                end
                
                % define modelID and categoryID
                curMID = ones(lenD_cur, 1) * j;
                curCID = ones(lenD_cur, 1) * i;
                
                category_id = [category_id; curCID];
                model_id = [model_id; curMID];
                
            end

        end
        i
    end


    lenF = length(list_el);
end

