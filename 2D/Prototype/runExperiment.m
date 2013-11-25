%> Name: runExperiment
%>
%> Description: The entry function to libHOP code. Calls libHOP with
%> desired parameters on specified dataset. For each image in the dataset,
%> we create a graph out of level 1 gabor filter responses, and form
%> edges. Then, we compress samples of each category with SUBDUE to get 
%>
%> @param datasetName Name of the dataset to work on. 
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 18.11.2013

function [] = runExperiment( datasetName )

    % Initial parameters
    numberOfFilters = 6;
    property = 'co-occurence';

    % Learn dataset path relative to this m file
    currentFileName = mfilename('fullpath');
    [currentFilePath, ~, ~] = fileparts(currentFileName);
    datasetFolder = [currentFilePath '/datasets/' datasetName '/'];
    addpath([currentFilePath '/utilities']);
    
    % Specify name of the graph file
    graphFileName = [currentFilePath '/graphs/' datasetName '.g'];
    fp = fopen(graphFileName, 'w');

    % Get all images under the dataset
    fileNames = fuf([datasetFolder '*_crop.png'], 1, 'detail');
    
    % Extract graphs of each image and add to the final graph
 %   for fileItr = 1:size(fileNames,1) 
    for fileItr = 1:20
        img = imread(fileNames{fileItr});
        nodes = getNodes(img, numberOfFilters, currentFilePath);
        edges = getEdges(nodes, property);
    
        % Print the graph to the file.
        printGraphToFile(fp, nodes, edges);
    end
    
    fclose(fp);
end

