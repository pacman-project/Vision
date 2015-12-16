folderName = 'D:\Input Data\VladislavSTD\Vladislav_STD\All reserves\depthAll\layer3';

file_list = load_filelist(folderName);
lenF = length(file_list);


for i = 2:lenF
    I = imread(file_list{i});
    imtool(I, [0, max(max(I))]);
    
end