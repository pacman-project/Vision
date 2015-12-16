% this script is to mearge two datasets
clear all;
dataSetMain = 'D:\PacmanVisionDataset\Views128\Dist_500\';
dataSetAdd = 'D:\PacmanVisionDatasetAdd\Views128\Dist_500\';

subfoldersAdd = list_of_subfolders(dataSetAdd);
len = length(subfoldersAdd);

for i = 3:len
    path = [dataSetMain, subfoldersAdd{i}];
    pathAdd = [dataSetAdd, subfoldersAdd{i}];
    subsubfolders = list_of_subfolders(path);
    cur = length(subsubfolders) - 1;
    
    subFoldersAdd =  list_of_subfolders(pathAdd);
    lenAdd = length(subFoldersAdd);
    for j = 3:lenAdd
        curFolderIn =  subFoldersAdd{j};
        pos = strfind(curFolderIn, '_');
        pos = pos(end);
        str = [curFolderIn(1:pos), num2str(cur)];
        strOut = [path,'\',str];
        strIn = [dataSetAdd, subfoldersAdd{i}, '\', curFolderIn];
        movefile(strIn, strOut);
        cur = cur + 1;
    end
end

a = 2;