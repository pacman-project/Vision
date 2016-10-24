%> Name: saveStats
%>
%> Description: Given the current graph and vocabulary, this function
%> calculates statistics such as shareability of nodes and image coverage
%> in terms of detected leaf nodes.
%>
%> @param vocabLevel Current vocabulary level.
%> @param graphLevel Current graph level.
%> @param leafNodes Array containing label, position and image id of leaf
%> nodes.
%> @param options Program options.
%> @param state Any of the following: 'preInhibition', 'postInhibition'.
%> The states are self-explanatory.
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 11.02.2014
function [avgShareability, avgCoverage, maxCoverageVals] = saveStats(vocabLevel, graphLevel, leafNodeCoords, maxCoverageVals, numberOfImages, options, state, levelItr)
    if isempty(vocabLevel) || isempty(graphLevel)
        avgShareability = 0;
        avgCoverage = 0;
        return;
    end
    
    % Find filter size, and area boundaries to cover. 
    filterSize = size(options.filters{1});
    filterSize = filterSize(1:2);
    topHalfSize = ceil(filterSize/2) - 1;
    bottomHalfSize = floor(filterSize/2);
    imageSize = options.imageSize;
    if levelItr == 1
         maxCoverageVals = zeros(numberOfImages,1);
    end
    leafNodeCoords = double(leafNodeCoords);
    
    %% Calculate shareability of each composition in vocabulary.
    % This is done by simply estimating number of images this composition is
    % seen in. 
    vocabLabelIds = [vocabLevel.label]';
    numberOfVocabNodes = max(vocabLabelIds);
    labelIds = vocabLabelIds([graphLevel.labelId]');
    imageIds = [graphLevel.imageId]';
    leafNodes = {graphLevel.leafNodes};
    compositionImageCountArr = zeros(numberOfVocabNodes,1);
    
    for labelId = 1:numberOfVocabNodes
       compositionImageCountArr(labelId) = numel(unique(imageIds(labelIds==labelId))) / numberOfImages;
    end
    avgShareability = mean(compositionImageCountArr(~isnan(compositionImageCountArr))) * 100;
    
    %% Calculate coverage of each image's leaf nodes in detected realizations.
    imageCoverageArr = zeros(numberOfImages,1);
    for imageItr = 1:numberOfImages
       catLeafNodes = [leafNodes{imageIds == imageItr}];
       if ~isempty(catLeafNodes)
           detectedLeafNodes = int32(fastsortedunique(sort([leafNodes{imageIds == imageItr}])));
       else
           continue;
       end
           
       detectedLeafNodeCoords = leafNodeCoords(detectedLeafNodes, :);
       imageMatrix = zeros(imageSize) > 0;
       for leafNodeItr = 1:numel(detectedLeafNodes)
            topLeft = detectedLeafNodeCoords(leafNodeItr,:) - topHalfSize;
            bottomRight = detectedLeafNodeCoords(leafNodeItr,:) + bottomHalfSize;
            imageMatrix(topLeft(1):bottomRight(1), topLeft(2):bottomRight(2)) = 1;
       end
       imageCoverageArr(imageItr) = nnz(imageMatrix);
       if levelItr == 1
           maxCoverageVals(imageItr) = nnz(imageMatrix);
       end
    end
    
    %% After obtaining coverage array, we calculate actual coverage in terms of image pixels.
    avgCoverage = mean(imageCoverageArr(imageCoverageArr > 0) ./ maxCoverageVals(imageCoverageArr > 0));
    
    %% Save info.
    statFile = [options.currentFolder '/output/' options.datasetName '/' state '_' num2str(levelItr) '.mat'];
    save(statFile, 'compositionImageCountArr', 'avgShareability', 'imageCoverageArr', 'avgCoverage');
    clearvars -except avgShareability avgCoverage maxCoverageVals
end

