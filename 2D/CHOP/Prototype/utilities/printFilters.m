%> Name: printFilters
%>
%> Description: Given square filters in a 
%>
%> @param W NxD array which contains N filters.
%> @param dimensions 1xM array, which contains the filter dimensions.
%> M is either 2 or 3.
%> @param datasetName The name of the dataset which this filters are
%> trained for. 
%>  
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 12.08.2015
function [] = printFilters( W, dimensions, datasetName )
    currentFileName = mfilename('fullpath');
    [currentPath, ~, ~] = fileparts(currentFileName);
    visFolderName = [currentPath '/filters/vis/' datasetName];
    filterFolderName = [currentPath '/filters/auto/' datasetName];
    
    save([options.currentFolder '/filters/vis/' datasetName '/C.mat'], 'C', 'mu', 'invMat', 'whMat');
    
end

