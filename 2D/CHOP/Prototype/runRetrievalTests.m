%> Name: runRetrievalTests
%>
%> Description: This function runs the main image retrieval tests to return
%> a confusion matrix, along with other performance values. What is 
%>
%> @param graphFolder The folder that includes each object's graphs. Each
%> object graph's path should be of the form:
%> graphFolder/objectClass/matchingLevel.g
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 20.12.2013
function [ ] = runRetrievalTests( currentFolder, datasetName, graphFolder, matchingLevel )
    objectGraphList = fuf([graphFolder '/*level' num2str(matchingLevel) '.g'], 1, 'detail');
    outputFolder = [currentFolder '/visualization/' datasetName];
    if ~exist(outputFolder, 'dir')
       mkdir(outputFolder); 
    end
%    n = numel(objectGraphList);
    costMatrix = zeros(numel(objectGraphList));
    for objectGraphItr = 1:numel(objectGraphList)
        for secGraphItr = 1:numel(objectGraphList)
            command = [currentFolder '/miners/gm ' objectGraphList{objectGraphItr} ' ' objectGraphList{secGraphItr}];
            [status, result] = system(command);
            % Parse result to get the match cost.
            if status == 0
                costIdx = strfind(result, 'Match Cost =');
                matchCost = textscan(result(1, (costIdx+13):end), '%f\n');
                matchCost = matchCost{1};
                costMatrix(objectGraphItr, secGraphItr) = matchCost;
            end
        end
    end
    matchMatrix = (costMatrix * -1) + max(max(costMatrix));
%    matchMatrix(1:(n+1):n*n) = 0;
    figure;
    colormap('hot');   % set colormap
    imagesc(matchMatrix);        % draw image and scale colormap to values range
    colorbar;          % show color scale
    saveas(gcf, [outputFolder '/matchMatrix.fig']);
    saveas(gcf, [outputFolder '/matchMatrix.png']);
    save([outputFolder '/matchMatrix.mat'], 'matchMatrix');
end

