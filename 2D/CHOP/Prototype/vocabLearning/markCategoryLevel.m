function [options, nodeCoverages, imageCoverages] = markCategoryLevel(graphLevel, level1Coords, levelItr, avgCoverage, options)
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
   if levelItr == options.categoryLevel || mean(imageCoverages(remainingImages)) >= options.categoryLevelCoverage
        options.stopFlag = true;
        options.categoryLevel = levelItr;
   end
end