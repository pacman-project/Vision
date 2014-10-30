function [finalImg] = printCloseFilters(filters, distanceMatrix)
    numberOfFilters = numel(filters);
    shownFilters = 6;
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
            filter1 = uint8(255 * (filter1 - min(min(min(filter1)))) / (max(max(max(filter1))) - min(min(min(filter1)))));
            finalImg(startX:(startX+imgDim(1)-1), startY:(startY+imgDim(1)-1), :) = filter1;
        end
    end
end