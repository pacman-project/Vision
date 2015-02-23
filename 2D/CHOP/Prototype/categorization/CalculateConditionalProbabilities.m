%> Name: CalculateConditionalProbabilities
%>
%> Description: this function works on the output of CHOP. It basically
%> lists the contributions of each node in the vocabulary to each category,
%> which helps with the categorization step, using the training images'
%> labels. The vocabulary is saved again to the same place after
%> processing.
%>
%> @param datasetName The name of the dataset.
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 03.07.2014
function [ ] = CalculateConditionalProbabilities( datasetName )
    %#ok<*NODEF>
    options = SetParameters(datasetName, 'train');
    % Read the vocabulary and the exported realizations. 
    load([options.currentFolder '/output/' datasetName '/vb.mat'], 'vocabulary', 'categoryNames');
    load([options.currentFolder '/output/' datasetName '/export.mat'], 'exportArr', 'categoryArrIdx', 'poseArr');
    
    % We go through each layer of the vocabulary, and update categoryArr of
    % every node in the vocabulary with the probabilities that this node
    % belongs to each category. The probabilities of each category for a
    % single node sum up to 1. categoryArr of a node is 1xN array, where N
    % is the number of categories.
    totalNumberOfParts = sum(cellfun(@(x) numel(x), vocabulary)); %#ok<USENS>
    probArr = zeros(totalNumberOfParts, 2 + numel(categoryNames));
    partOffset = 1;
    for levelItr = 1:numel(vocabulary)
        
        vocabLevel = vocabulary{levelItr};
        
        % Put exported realizations in different cells in a cell array 
        % (for fast parallel processing).
        realizations = exportArr(exportArr(:,4) == levelItr, :); 
        numberOfParts = numel(vocabLevel);
        realArr = cell(numberOfParts,1);
        probArr(partOffset:(partOffset+numberOfParts-1),1) = levelItr;
        probArr(partOffset:(partOffset+numberOfParts-1),2) = 1:numberOfParts;
        
        for partItr = 1:numberOfParts
            realArr{partItr} = categoryArrIdx(unique(realizations(realizations(:,1) == partItr, 5)));
            for categoryItr = 1:numel(categoryNames)
               probArr(partOffset + partItr - 1, 2 + categoryItr) = ...
                   nnz(realArr{partItr} == categoryItr)/numel(categoryArrIdx == categoryItr) * 100;
            end
        end
        
        
        % Save probabilities.
        partOffset = partOffset + numel(vocabLevel);
    end
    
    % Save vocabulary.
    save([options.currentFolder '/output/' datasetName '/prob.mat'], 'probArr');
end

