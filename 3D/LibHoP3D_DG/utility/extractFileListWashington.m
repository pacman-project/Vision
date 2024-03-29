% this is to extract filelist fromt the Washington dataset

function [list_depth, list_mask, list_RGB, lenF] = extractFileListWashington(depthPath, is_subset, subsetPercent)

    disp('extracting file list...');
    list_depth = {};
    list_mask = {};
    list_RGB = {};

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
            list_file = load_filelist(strS);
            list_file = list_file';
            lenSSF = length(list_file);

            indsD = 2:3:lenSSF;
            indsRGB = 1:3:lenSSF;
            indsMask = 3:3:lenSSF;

            fileD = list_file(indsD);
            fileRGB = list_file(indsRGB);
            fileMask = list_file(indsMask);

            lenSSF3 = length(fileD); % in fact lenSSF/3

            % select images randomly

            if (is_subset == true)
                if (subsetPercent > 1) || (subsetPercent < 0)
                    disp('Error: subset percentage is invalid');
                    %return;
                end
                subset_len = round(lenSSF3 * subsetPercent);
                % select the subset
                nums = randperm(lenSSF3);
                nums = nums(1:subset_len);

                list_depth = [list_depth; fileD(nums)];
                list_mask =  [list_mask;  fileMask(nums)];
                list_RGB =   [list_RGB;   fileRGB(nums)];
            else
                list_depth = [list_depth; fileD];
                list_mask =  [list_mask;  fileMask];
                list_RGB =   [list_RGB;   fileRGB];
            end

        end
        i

    end


    lenF = length(list_depth);
end

