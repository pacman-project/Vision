%> Name: assignEdgeLabels
%>
%> Description: Given modes, and edges in graphLevel, this function assigns
%> an updated edge label to all edges in graphLevel, based on the learned
%> modes of the 2D distributions in 'modes' argument.
%> 
%> @param vocabLevel Vocabulary level.
%> @param graphLevel The main graph that includes the nodes and their
%> respective edges.
%> @param modes The mode list representing edge categories.
%>               modes are of the form: [ nodeLabel1, nodeLabel2, edgeId,
%> rfCoord1, rfCoord2, cov11, cov12, cov21, cov22, coord1, coord2, weight;
%>  ...];
%> @param edgeCoords The coordinates of all possible edge types.
%> 
%> @retval graphLevel The main graph with updated edge labels.
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 06.09.2015
function [ graphLevel ] = assignEdgeLabels(graphLevel, modes, modeProbArr, edgeCoords)
    display('Assigning edge labels using modes...');
     allEdges = cat(1, graphLevel.adjInfo);
     nodeIds = [graphLevel.labelId]';
     halfSize = ceil(size(modeProbArr,2)/2);
     edgeCoords = edgeCoords + halfSize;
     
     % If there are no edges, return.
     if isempty(allEdges)
          return;
     end
     
     [~, edgeOffsets] = ismember(1:numel(graphLevel), allEdges(:,1));
     edgeOffsets = edgeOffsets - 1;
     
     % Create mode index array.
     [uniqueModes, IA, ~] = unique(modes(:,1:2), 'stable', 'rows');
     maxMode = numel(IA);
     uniqueModes = int32(uniqueModes);
     
     [~, modeMapIds] = ismember(nodeIds(allEdges(:, 1:2)), uniqueModes, 'rows');
     allEdges = {graphLevel.adjInfo};
     
%     totalEdgeCount = 0;
%     validEdgeCount = 0;
     
     % Get node ids of edges.
     zeroEdges = zeros(0, 'int32');
     modes = modes(:,3);
     assignedEdges = cell(numel(graphLevel),1);
     parfor graphLevelItr = 1:numel(graphLevel)
          edges = allEdges{graphLevelItr};

          % Edges empty, do nothing.
          if isempty(edges)
               assignedEdges{graphLevelItr} = zeroEdges;
               continue;
          end
          
          edgeOffset = edgeOffsets(graphLevelItr);
          
          % Assign each edge a to a relevant mode.
          numberOfEdges = size(edges,1);
          for edgeItr = 1:numberOfEdges
               % Get the mode map id.
               modeRow = modeMapIds(edgeOffset + edgeItr);
               
               if modeRow == 0
                  continue; 
               end
               
               % Get possible row ids.
               if modeRow == maxMode
                   relevantModes = modes(IA(modeRow):end, 1);
               else
                   relevantModes = modes(IA(modeRow):(IA(modeRow+1)-1), 1);
               end
               
               % Finally, obtain the edge id.
               relevantCoords = edgeCoords(edges(edgeItr,3),:);
               clusterId = modeProbArr(relevantCoords(1),relevantCoords(2), modeRow);
               if clusterId == 0
                    newLabel = 0;
               else
                    newLabel = int32(relevantModes(clusterId));
               end
               
               % Assign new label.
               edges(edgeItr,3) = newLabel;
          end
          % Collect statistics about what percentage of edges are rendered
          % useless.
%          totalEdgeCount = totalEdgeCount + size(edges,1);
%          validEdgeCount = validEdgeCount + nnz(edges(:,3) > 0);
          
          % Save edges.
          assignedEdges{graphLevelItr} = edges(edges(:,3) > 0, :);
     end
     
     % Assign the edges to the relevant data structure.
     [graphLevel.adjInfo] = deal(assignedEdges{:});
     
     % Report on deleted edges.
%     display(['%' num2str(round(100*(validEdgeCount/totalEdgeCount)))...
%         ' of all edges have been preserved. Total: ' num2str(totalEdgeCount) ...
%         ' edges. Deleted:' num2str(totalEdgeCount - validEdgeCount) ' edges.']); 
     
     clearvars -except graphLevel
end