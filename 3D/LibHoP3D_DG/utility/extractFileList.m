% this is to extract filelist and subset of it (if required)

function [list_depth, len] = extractFileList(depthPath, is_subset, subsetPercent)

    list_depth = load_filelist(depthPath);

    len = length(list_depth);
    subset_len = round(len * subsetPercent);

    if (is_subset == true)
        if (subset_len > len)
            disp('Error: subset_len is larger than list of files');
            return;
        end
        % select the subset
        nums = randperm(len);
        nums = nums(1:subset_len);
        list_depth = list_depth(nums);
        len = subset_len;
    end


end

