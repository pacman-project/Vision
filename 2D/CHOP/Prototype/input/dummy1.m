% files = fuf([pwd '/input/ethz/*.pgm'], 1, 'detail');
% for fileItr = 1:numel(files)
%    img = imread(files{fileItr});
%    [filePath, fileName, ~] = fileparts(files{fileItr});
% %   img = img>0;
%    imwrite(img, [filePath '/' fileName '.png']);
% end
% images=loadMNISTImages('./input/vocab/MNIST/train-images.idx3-ubyte')';
% labels = loadMNISTLabels('./input/vocab/MNIST/train-labels.idx1-ubyte');
% uniqueLabels = unique(labels);
% for labelItr = 1:numel(uniqueLabels)
%     if ~exist(['./input/vocab/MNIST/' num2str(uniqueLabels(labelItr))], 'file')
%         mkdir(['./input/vocab/MNIST/' num2str(uniqueLabels(labelItr))]);
%     end
% end
% for itr = 1:size(images,1);
%     image = reshape(images(itr,:), [28 28]);
%     imwrite(image, ['./input/vocab/MNIST/' num2str(labels(itr)) '/' num2str(itr) '.png']);
% endi
images=loadMNISTImages('./input/MNIST/t10k-images.idx3-ubyte')';
labels = loadMNISTLabels('./input/MNIST/t10k-labels.idx1-ubyte');
uniqueLabels = unique(labels);
for labelItr = 1:numel(uniqueLabels)
    if ~exist(['./input/MNIST/test/' num2str(uniqueLabels(labelItr))], 'file')
        mkdir(['./input/MNIST/test/' num2str(uniqueLabels(labelItr))]);
    end
end
for itr = 1:size(images,1);
    image = reshape(images(itr,:), [28 28]);
    imwrite(image, ['./input/MNIST/test/' num2str(labels(itr)) '/' num2str(itr) '.png']);
end