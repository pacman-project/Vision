% extract folder name given full path

function [folderName, fileName] = getFolderName(path)
    pos1 = strfind(path, '/');
    pos2 = strfind(path, '\');
    pos1 = [pos1, pos2];
    mp = max(pos1);
    folderName = path(1:mp);
    fileName = path(mp+1:end);
end