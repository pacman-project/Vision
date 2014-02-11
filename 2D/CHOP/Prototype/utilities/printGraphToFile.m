%> Name: printGraphToFile
%>
%> Description: Given a set of nodes and edges, we form the graph and print
%> it to the given file.
%>
%> @param fp File handle
%> @param nodes Nodes of the graph
%> @param edges Edges of the graph
%> @param isNumericNode 1 if node labels are numeric, 0 otww
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 20.11.2013
function [ ] = printGraphToFile( fp, nodes, edges, isNumericNode )
    
    % Print nodes
    if isNumericNode
        for nodeItr = 1:size(nodes,1)
            fprintf(fp, ['v ' num2str(nodeItr) ' ' num2str(nodes(nodeItr)) '\n']);
        end
    else
        for nodeItr = 1:size(nodes,1)
            fprintf(fp, ['v ' num2str(nodeItr) ' ' nodes(nodeItr) '\n']);
        end
    end
    
    % Print edges
    for edgeItr = 1:size(edges,1)
        if edges(edgeItr,4) == 1
            edgeId = 'd';
        else
            edgeId = 'u';
        end
        fprintf(fp, [edgeId ' ' num2str(edges(edgeItr,1)) ' ' num2str(edges(edgeItr,2)) ' ' num2str(edges(edgeItr,3)) '\n']); 
    end
end

