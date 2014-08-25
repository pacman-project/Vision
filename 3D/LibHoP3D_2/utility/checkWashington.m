function [i,j, is_end] = checkWashington(startI, startK)

    % this is a script to check the Washington dataset  

    depthPath = '/home/vvk201/Wash-rgbd-dataset';

    is_end = false;
    %  extract a list of subfolders
    subfolders = list_of_subfolders(depthPath);
    % first two elements are garbage
    subfolders(1:2) = [];
    lenS = length(subfolders);

    for i = startI:lenS  % for each category
        str = [depthPath, '/', subfolders{i}];
        % select a list of subfolders again
        subSubFolders = list_of_subfolders(str);
        subSubFolders(1:2) = [];
        lenSS = length(subSubFolders);

        for j = startK:lenSS  % for each object

            strS = [str, '/', subSubFolders{j}];
            list_file = load_filelist(strS);
            list_file = list_file';
            lenSSF = length(list_file);

            is_successfull = false;

            for k = 1:3:lenSSF
                imageName = list_file{k};
                imageDepth =  list_file{k+1};
                imageMask =  list_file{k+2};

                imageName = imageName(1:end - 4 - 4);
                imageDepth = imageDepth(1:end - 9 - 4);
                imageMask = imageMask(1:end - 8 - 4);

                if length(imageName) ~= length(imageDepth) || length(imageName) ~= length(imageMask)
                    disp(imageName);
                    disp(imageDepth);
                    disp(imageMask);
                    filename = [imageName, 'loc.txt'];
                    delete(list_file{k}, list_file{k+1}, filename);
                    return;
                end

                b1 = max(imageDepth ~= imageMask);
                b2 = max(imageName ~= imageDepth);

    %             if (b1 + b2 > 0)
    %                 str = [imageName, ' ', imageDepth, ' ', imageDepth];
    %                 disp(str);
    %                 a = 2;
    %             end
            end
        end
        i
    end
    is_end = true;
end

