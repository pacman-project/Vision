%> Name: printVocabulary
%>
%> Description: Given the vocabulary and low-level filters, this function
%> visualizes each node in the vocabulary across all layers based on how they
%> are represented. 
%>
%> @param datasetName The name of the dataset for which we'll visualize the
%> vocabulary.
%> 
%> @retval none
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 26.02.2015
function [ ] = printVocabulary( datasetName )
    % Get program options and relevant parameters.
    options = SetParameters(datasetName, false);
    outputFolder = options.outputFolder;
    debugFolder = options.debugFolder;
    threshold = options.subdue.threshold;
    newVisFolder = [debugFolder '/vocabVisualization'];
    if ~exist(newVisFolder, 'dir')
       mkdir(newVisFolder); 
    end
    
    % Read vocabulary.
    load([outputFolder '/vb.mat'], 'vocabulary');
    filters = options.filters;
     
    %% Starting from bottom-up, we will print vocabulary nodes one by one.
    for levelItr = 1:numel(vocabulary)
       vocabLevel = vocabulary{levelItr}; 
       
       % If first level, we'll simply print the filters.
       if levelItr == 1
           
           
       end
        
    end
end

