%> Name: visualizeLevel
%>
%> Description: Given the current and previous vocabulary level, this
%> function visualizes the current vocabulary level as a separate patche for
%> each word in the vocabulary level.
%> 
%> @param currentLevel Current vocabulary level.
%> @param graphLevel Current graph level.
%> @param levelId Identifier of the current level.
%> @param numberOfPrevNodes Number of nodes in previous vocabulary level.
%> @param options Program options.
%> @param isRedundant If currentLevel consists of redundant compositions,
%> set 1. Otherwise set 0.
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 10.02.2014
%> Redundant vocabulary output option added. 10.05.2014
function [] = visualizeORNodes( currentLevel, levelId, options)
    % Read options to use in this file.
    currentFolder = options.currentFolder;
    datasetName = options.datasetName;
    labelIds = double(cat(1, currentLevel.label));
    reconstructionDir = [options.debugFolder '/level' num2str(levelId) '/modalProjection/'];
    load([reconstructionDir 'muImgs.mat']);
    maskSize = size(smallMuImgs,2);
    
    % Get the count of the most frequent label.
    if numel(unique(labelIds)) > 1
         labelCounts = hist(labelIds, double(unique(labelIds)));
         maxCount = max(labelCounts);
         [~, sortIdx] = sort(labelCounts, 'descend');
    else
         maxCount = numel(labelIds);
         sortIdx = 1:numel(labelIds);
    end
    
    % Using the maximum dimensions, transform each composition image to the
    % same size. ß
    overallImage = ones(numel(unique(labelIds)) * (maskSize+1)+1, (maxCount)*(maskSize+1)+1, 1, 'uint8') * 255;

    % Write images one by one.
    for labelItr = 1:numel(sortIdx)
         optionSet = find(labelIds == sortIdx(labelItr));
         rowStart = 2 + (labelItr-1)*(maskSize+1);
         for optionItr = 1:numel(optionSet)
              colStart = 2 + (optionItr-1) * (maskSize+1);
              sampleImg = squeeze(smallMuImgs(optionSet(optionItr), :, :));
              overallImage(rowStart:(rowStart+maskSize-1), ...
                    colStart:(colStart+maskSize-1), :) = sampleImg;
         end
    end
    
    % A final make up in order to separate masks from each other by 1s.
    imwrite(overallImage, [currentFolder '/debug/' datasetName '/level' num2str(levelId) '_ORNodeSets.jpg']);
    imwrite(overallImage, [currentFolder '/debug/' datasetName '/level' num2str(levelId) '_ORNodeSets_HQ.png']);
    clearvars
end



