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
function [avgShareability, avgCoverage] = saveStats(vocabLevel, graphLevel, leafNodes, numberOfImages, options, state, levelItr)
    %% Calculate shareability of each composition in vocabulary.
    % This is done by simply assessing number of images this composition is
    % seen in. 
    labelIds = [graphLevel.labelId]';
    imageIds = [graphLevel.imageId]';
    compositionImageCountArr = zeros(numel(vocabLevel),1);
    
    for labelId = 1:numel(vocabLevel)
       compositionImageCountArr(labelId) = numel(unique(imageIds(labelIds==labelId))) / numberOfImages;
    end
    avgShareability = mean(compositionImageCountArr);
    
    %% Calculate coverage of each image's leaf nodes in detected realizations.
    detectedLeafNodes = unique([graphLevel.leafNodes]);
    leafNodeIds = 1:size(leafNodes,1);
    leafImageIds = cell2mat(leafNodes(:,3));
    overallCoverage = numel(detectedLeafNodes) / size(leafNodes,1);
    imageCoverageArr = zeros(numberOfImages,1);
    for imageItr = 1:numberOfImages
        imageCoverageArr(imageItr) = numel(intersect(detectedLeafNodes, leafNodeIds(leafImageIds == imageItr))) / ...
            numel(leafNodeIds(leafImageIds == imageItr));
    end
    avgCoverage = mean(imageCoverageArr);
    
    %% Save info.
    statFile = [options.currentFolder '/output/' options.datasetName '/' state '_' num2str(levelItr) '.mat'];
    save(statFile, 'compositionImageCountArr', 'avgShareability', 'overallCoverage', 'imageCoverageArr', 'avgCoverage');
end

