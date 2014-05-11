function [nameFolds] =  list_of_subfolders(pathFolder)

    d = dir(pathFolder);
    isub = [d(:).isdir]; %# returns logical vector
    nameFolds = {d(isub).name}';
end