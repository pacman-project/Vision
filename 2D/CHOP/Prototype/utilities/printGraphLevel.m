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
    
    %% For each image, select the nodes and their edges, and print them.
    nodeOffset = 0;
    rfIds = cat(1, graphLevel.rfId);
    uniqueRfIds = unique(rfIds)';
    
    for imageId = uniqueImageIds
        imageIdx = imageIds==imageId;
        for rfId = uniqueRfIds
            rfIdx = rfIds==rfId & imageIdx;
            if nnz(rfIdx) < 1
                continue;
            end
            
            rfNodeIds = cat(1, graphLevel(rfIdx).labelId);
            rfEdges = cat(1, graphLevel(rfIdx).adjInfo);
            
            if numel(rfNodeIds) < 1
                continue;
            end
            if nnz(rfEdges) > 0
                rfEdges(:,1:2) = rfEdges(:,1:2) - nodeOffset;
            end

            % Print the graph to the file.    
            % Print positive graph indicator
            fprintf(fp, 'XP\n');
            printGraphToFile(fp, rfNodeIds, rfEdges, true);
            
            nodeOffset = nodeOffset + numel(rfNodeIds);
        end
    end
    fclose(fp);
end
    
        