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
    sampleImg = imread([reconstructionDir num2str(1) '.png']);
    compMaskSize = size(sampleImg);
    
    % Get the count of the most frequent label.
    if numel(unique(labelIds)) > 1
         labelCounts = hist(labelIds, double(unique(labelIds)));
         maxCount = max(labelCounts);
         [~, sortIdx] = sort(labelCounts, 'descend');
    else
         labelCounts = numel(labelIds);
         maxCount = numel(labelIds);
         sortIdx = 1:numel(labelIds);
    end
    
    % Using the maximum dimensions, transform each composition image to the
    % same size. 
    if numel(compMaskSize) == 2
         bandCount = 1;
    else
         bandCount = compMaskSize(3);
    end
    overallImage = NaN(numel(unique(labelIds)) * (compMaskSize(2)+1)+1, (maxCount)*(compMaskSize(1)+1)+1, bandCount);

    % Write images one by one.
    for labelItr = 1:numel(sortIdx)
         optionSet = find(labelIds == sortIdx(labelItr));
         rowStart = 2 + (labelItr-1)*(compMaskSize(1)+1);
         for optionItr = 1:numel(optionSet)
              colStart = 2 + (optionItr-1) * (compMaskSize(2)+1);
              sampleImg = imread([reconstructionDir num2str(optionSet(optionItr)) '.png']);
              overallImage(rowStart:(rowStart+compMaskSize(1)-1), ...
                    colStart:(colStart+compMaskSize(2)-1), :) = sampleImg;
         end
    end
    
    % A final make up in order to separate masks from each other by 1s.
    overallImage(isnan(overallImage)) = 0;
    overallImage = uint8(overallImage);

    % Then, write the compositions the final image.
    imwrite(overallImage, [currentFolder '/debug/' datasetName '/level' num2str(levelId) '_ORNodeSets.png']);
    clearvars
end



