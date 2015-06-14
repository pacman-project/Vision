folderName = 'F:\Princeton Rendered\All Views';
list_depth = load_filelist(folderName);

lenFolderName = length(folderName);
outFolder = 'F:\Princeton Rendered\All Views2';

lenF = length(list_depth);

parfor i = 1:lenF
    
    str = list_depth{i};
    I = imread(str);
    
    filename = str(lenFolderName+1:end);
    filename(end-4) = '2';
    outName = [outFolder, filename];
    
    %I = flipud(I);
    I = imrotate(I, 90);
    imwrite(I, outName, 'png');
    
    if mod(i, 100) == 0
        i
    end
end