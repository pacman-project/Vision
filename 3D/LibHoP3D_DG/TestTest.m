

depthPath = 'D:\Input Data\VladislavSTD\Vladislav_STD\depthEdgesCorners_layer2';

list_depth = load_filelist(depthPath);

for i = 1:length(list_depth);
     I = imread(list_depth{i});
     imtool(I, [0,1]);
     
     a = 2;
end