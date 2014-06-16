fileNames = fuf([pwd '/input/Swans/*.pgm'], 1, 'detail');
for fileItr = 1:numel(fileNames)
    [filePath, fileName, ext] = fileparts(fileNames{fileItr});
    tempImg = imread(fileNames{fileItr});
    tempImg = tempImg ~= median(median(double(tempImg)));
    imwrite(tempImg, [filePath '/' fileName, '.png']);
end
