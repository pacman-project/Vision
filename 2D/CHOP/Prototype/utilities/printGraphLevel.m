%> Name: printGraphLevel
%>
%> Description: Print the graph level to the file, as a set of nodes and
%their corresponding edges.
%> 
%> @param graphFileName The file to print in.
%> @param graphLevel The object graphs' current level to be printed.
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 03.02.2014
function [] = printGraphLevel(graphFileName, graphLevel)
    %% Get relevant info from graph level.
    imageIds = cat(1, graphLevel.imageId);
    uniqueImageIds = unique(imageIds)';
    fp = fopen(graphFileName, 'w');
    
    nodeSets = cell(numel(uniqueImageIds),1);
    for imageItr = 1:numel(uniqueImageIds)
        nodeSets(imageItr) = {graphLevel(imageIds==uniqueImageIds(imageItr))};
    end
    
    %% For each image, select the nodes and their edges, and print them.
    for imageItr = 1:numel(uniqueImageIds)
        nodeOffset = numel(find(imageIds<uniqueImageIds(imageItr)));
        imageNodes = nodeSets{imageItr};
        if isempty(imageNodes)
            continue;
        end

        imageNodeIds = cat(1, imageNodes.labelId);
        imageEdges = cat(1, imageNodes.adjInfo);

        if nnz(imageEdges) > 0
            imageEdges(:,1:2) = imageEdges(:,1:2) - nodeOffset;
        end

        % Print the graph to the file.    
        % Print positive graph indicator
        fprintf(fp, 'XP\n');
        printGraphToFile(fp, imageNodeIds, imageEdges, true);

    end
    fclose(fp);
end
    
        