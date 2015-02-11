% before running this one run those two scripts:
% addPartIdToDetections -- to obtain partIds
% addCategoryToDetections -- to obtain categoryIds
% these scripts are time consuming so you can load precomputed data from
% mat files: partIds.mat and categoryIds.mat
numCategories = max(categoryIds);
numParts = max(partIds);
partFreq = zeros(numParts, numCategories);
for i = 1:numParts
    disp(i);
        inds = find(partIds == i);
        for j=1:length(inds)
            partFreq(i,categoryIds(inds(j))) = partFreq(i,categoryIds(inds(j)))+1;
        end
end