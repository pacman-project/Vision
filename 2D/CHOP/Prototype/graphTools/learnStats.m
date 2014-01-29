%> Name: learnStats
%>
%> Description: Given the object graphs in training set, this function
%> automatically learns the geometric statistics between nodes and keeps them in
%> 'modes' matrix. High-level relations are also learned, and they are
%> returned in 'highLevelModes' array.
%> 
%> @param nodes Nodes from the current level.
%> @param mainGraph The object graphs' data structure.
%> @param leafNodeAdjArr Leaf node adjacency array, which is basically an
%>      adjacency list of every leaf node that is connected because they are in
%>      the neighborhood of each other.
%> @param options Program options.
%> @param currentLevel The current scene graph level number. Used to index
%>      right object graphs from mainGraph.
%> @param datasetName Name of the dataset.
%> 
%> @retval currentModes The mode list representing edge categories for this level. 
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 28.01.2014
function [currentModes, currentHighLevelModes] = learnStats(nodes, mainGraph, leafNodeAdjArr, options, currentLevel, datasetName)
    %% Step 1: Learn geometric relations between nodes belonging to a single object.
    if strcmp(options.property, 'mode')
        %% Learn modes from the data.
        [currentModes] = learnModes(nodes, mainGraph, leafNodeAdjArr, options, currentLevel, datasetName); 
    elseif strcmp(options.property, 'hist')
        %% In 'hist' type geometric property, we assign region centers in hMatrix as modes.
        % Please note that the output modes do not include any composition
        % id, and are therefore valid for all possible pairs.
        scale = (1/options.scaling)^(currentLevel-1);
        neighborhood = fix(options.edgeRadius * scale);
        
        % Read centers of regions in the histogram.
        load('hMatrix.mat', 'hMatrix'); 
        centroids = regionprops(hMatrix', 'Centroid');
        halfSizeMatrix = round(size(hMatrix,1)/2);
        centers = [centroids.Centroid];
        centers = [centers(1:2:end)', centers(2:2:end)'];
        
        % Normalize the center positions by translating the centers such
        % that midpoint of the matrix is the origin, and the they are
        % resized to fit current neighborhood.
        normCenters = round((centers - halfSizeMatrix) * (neighborhood/halfSizeMatrix));
        
        % Keep centers of histogram 
        currentModes = [zeros([size(normCenters,1), 2]), normCenters];
    else
        %% In co-occurence type property, there is no mode info.
        currentModes = [];
    end
    
    %% Step 2: Learn high level relations between nodes belonging to DIFFERENT objects.
    % These relations may correspond to high-level edges representing
    % view-related or category-related relations.
    % TODO: Explain format of these modes here.
    if strcmp(options.highLevelProperty, 'multiview')
        % TODO: Multi-view related modes will be formed.
        currentHighLevelModes = [];
    elseif strcmp(options.highLevelProperty, 'category')
        % TODO: Add category-wise modes. 
        currentHighLevelModes = [];
    else
        % 'none' option: Do not learn any high-level modes here.
        currentHighLevelModes = [];
    end 
end