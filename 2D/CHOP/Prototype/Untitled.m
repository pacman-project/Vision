renderFileNames = fuf([pwd '/render/*.png'], 1, 'detail');
for renItr = 1:numel(renderFileNames)
    renImg = imread(renderFileNames{renItr});
    renImg = imresize(renImg, 0.75);
%    renImg2 = rgb2gray(renImg)>0;
%    ind = find(renImg2);
%    [xInd, yInd] = ind2sub(size(renImg2), ind);
%    renImg = renImg(min(xInd):max(xInd), min(yInd):max(yInd), :);
   imwrite(renImg, renderFileNames{renItr});
end