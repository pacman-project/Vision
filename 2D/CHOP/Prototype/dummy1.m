files = fuf([pwd '/input/ethz/*.pgm'], 1, 'detail');
for fileItr = 1:numel(files)
   img = imread(files{fileItr});
   [filePath, fileName, ~] = fileparts(files{fileItr});
%   img = img>0;
   imwrite(img, [filePath '/' fileName '.png']);
end