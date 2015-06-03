% this is to extract filelist fromt the Washington dataset

function [list_depth, list_mask, list_RGB, lenF] = extractFileListMirela(depthPath, is_subset, subsetPercent, is_RGB, is_Mask)
    
    if nargin == 3
        is_RGB = true;
        is_Mask = false;
    end

    disp('extracting file list...');
    list_depth = {};
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

    for i = 1:lenS % for each category
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

                list_depth = [list_depth; fileD(nums)];
                
                if is_RGB
                    list_RGB =   [list_RGB;   fileRGB(nums)];
                end
                if is_Mask
%                   list_mask =  [list_mask;  fileMask(nums)];
                end
                
            else
                list_depth = [list_depth; fileD];
                if is_RGB
                    list_RGB =   [list_RGB;   fileRGB];
                end
                if is_Mask
%                 list_mask =  [list_mask;  fileMask];
                end
                
            end

        end
        i
    end


    lenF = length(list_depth);
end

