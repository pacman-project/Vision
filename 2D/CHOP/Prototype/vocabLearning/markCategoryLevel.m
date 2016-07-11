function [options, nodeCoverages, imageCoverages] = markCategoryLevel(graphLevel, level1Coords, avgCoverage, levelItr, options)
   % Allocate space for data structures.
   imageIds = [graphLevel.imageId];
   leafNodeCounts = cellfun(@(x) numel(x), {graphLevel.leafNodes});
   remainingImages = fastsortedunique(sort(imageIds));
   maxRemImg = double(max(imageIds));
   imageRemainingLeafCounts = hist(level1Coords(:,1), 1:maxRemImg);
   imageCoverages = zeros(maxRemImg, 1);
   nodeCoverages = (leafNodeCounts ./ imageRemainingLeafCounts(imageIds))/avgCoverage;
   
   % Find max coverage for every image.
   parfor imageItr = 1:maxRemImg
       idx = imageIds == imageItr;
       if nnz(idx) > 0
           imageCoverages(imageItr) = max(nodeCoverages(idx));
       end
   end

   % If conditions are met, move on.
   if nnz(imageCoverages(remainingImages) >= options.categoryLevelCoverage) == nnz(imageCoverages) || levelItr == options.categoryLevel
        options.stopFlag = true;
   end
end