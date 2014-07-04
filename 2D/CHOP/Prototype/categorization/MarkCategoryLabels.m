%> Name: MarkCategoryLabels
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
function [ ] = MarkCategoryLabels( datasetName )
    %#ok<*NODEF>
    options = SetParameters(datasetName, 'train');
    % Read the vocabulary and the exported realizations. 
    load([options.currentFolder '/output/' datasetName '/vb.mat'], 'vocabulary');
    load([options.currentFolder '/output/' datasetName '/export.mat'], 'exportArr', 'categoryArrIdx', 'categoryNames', 'poseArr');
    
    % We go through each layer of the vocabulary, and update categoryArr of
    % every node in the vocabulary with the probabilities that this node
    % belongs to each category. The probabilities of each category for a
    % single node sum up to 1. categoryArr of a node is 1xN array, where N
    % is the number of categories.
    for levelItr = 1:numel(vocabulary)
        vocabLevel = vocabulary{levelItr};
        
        % Put exported realizations in different cells in a cell array 
        % (for fast parallel processing).
        realizations = exportArr(exportArr(:,4) == levelItr, :); 
        numberOfParts = max(realizations(:,1));
        numberOfCategories = numel(categoryNames);
        realArr = cell(numberOfParts,1);
        
        for partItr = 1:numberOfParts
            realArr{partItr} = categoryArrIdx(realizations(realizations(:,1) == partItr, 5));
        end
        
        % Process each node in the vocabulary.
        probArr = cell(numberOfParts,1);
        for partItr = 1:numel(vocabLevel)
            probArr{partItr} = hist(realArr{partItr}, 1:numberOfCategories) / numel(realArr{partItr});
        end
        
        % Save probabilities.
        [vocabLevel.categoryArr] = deal(probArr{:});
        vocabulary{levelItr} = vocabLevel;
    end
    
    % Save vocabulary.
    save([options.currentFolder '/output/' datasetName '/vb.mat'], 'vocabulary');
end

