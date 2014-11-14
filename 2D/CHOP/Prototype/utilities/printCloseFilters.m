function [finalImg] = printCloseFilters(distanceMatrix, levelItr, options)
    if levelItr == 1
        filters = options.filters;
    else
        filters = cell(size(distanceMatrix,1),1);
        for filterItr = 1:numel(filters)
            datasetName = options.datasetName;
            filters{filterItr} = imread([options.currentFolder '/debug/' datasetName '/level' num2str(levelItr) '/reconstruction/' num2str(filterItr) '_uni.png']);
        end
    end
    for filterItr = 1:numel(filters)
       filePath = [options.currentFolder '/debug/' options.datasetName '/level' num2str(levelItr) '/cropped/' num2str(filterItr) '_1.png'];
       if exist(filePath, 'file')
           orgImg = imread(filePath);
       else
           orgImg = [];
       end
       filter1 = double(filters{filterItr});
       filter1 = uint8(255 * (filter1 - min(min(min(filter1)))) / (max(max(max(filter1))) - min(min(min(filter1)))));
%        if ~isempty(orgImg) && size(filter1,1) == size(orgImg,1) && size(filter1,2) == size(orgImg,2)
%            newFilterImg = zeros(size(filter1,1), size(filter1,2), size(orgImg,3), 'uint8');
%            if size(filter1,3) < size(orgImg,3)
%                for bandItr = 1:size(orgImg,3)
%                     newFilterImg(:,:,bandItr) = uint8(round(double(filter1(:,:,1)) * 0.7 + double(orgImg(:,:,bandItr)) * 0.3));
%                end
%            else
%                for bandItr = 1:size(orgImg,3)
%                     newFilterImg(:,:,bandItr) = uint8(round(double(filter1(:,:,bandItr)) * 0.7 + double(orgImg(:,:,bandItr)) * 0.3));
%                end
%            end
%            filter1 = newFilterImg;
%        end
       filters{filterItr} = filter1;
    end
    numberOfFilters = numel(filters);
    shownFilters = min(size(distanceMatrix,1), 6);
    imgDim = size(filters{1});
    finalImg = zeros([numberOfFilters + 1 + numberOfFilters*imgDim(1), shownFilters + 1 + imgDim(1) * shownFilters, size(filters{1},3)], 'uint8');
    for filterItr = 1:numel(filters)
        scores = distanceMatrix(filterItr,:);
        [~, rankings] = sort(scores, 'ascend');
        rankings = rankings(1:shownFilters);
        for shownItr = 1:shownFilters
            startX = 2 + (filterItr-1) * (imgDim(1)+1);
            startY = 2 + (shownItr-1) * (imgDim(1)+1);
            filter1 = filters{rankings(shownItr)};
            finalImg(startX:(startX+imgDim(1)-1), startY:(startY+imgDim(2)-1), :) = filter1;
        end
    end
    imwrite(finalImg, [options.currentFolder '/debug/' options.datasetName '/level' num2str(levelItr) '_neighbors.png']);
end