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
function [ graphLevel ] = assignEdgeLabels(vocabLevel, graphLevel, modes, edgeCoords)
    display('Assigning edge labels using modes...');
     allEdges = cat(1, graphLevel.adjInfo);
     vocabLevelLabels = [vocabLevel.label]';
     nodeIds = vocabLevelLabels([graphLevel.labelId]');
     
     % If there are no edges, return.
     if isempty(allEdges)
          return;
     end
     
     % Get node ids of edges.
     parfor graphLevelItr = 1:numel(graphLevel)
          edges = graphLevel(graphLevelItr).adjInfo;

          % Edges empty, do nothing.
          if isempty(edges)
               continue;
          end
          
          edgeNodeLabels = nodeIds(edges(:, 1:2));
          if size(edgeNodeLabels,2) ~= 2
               edgeNodeLabels = edgeNodeLabels';
          end
          
          % Assign each edge a to a relevant mode.
          validEdgeIdx = ones(size(edges,1),1) > 0;
          for edgeItr = 1:size(edges,1)
               
               relevantModes = modes(modes(:, 1) == edgeNodeLabels(edgeItr,1) & modes(:,2) == edgeNodeLabels(edgeItr,2), :);
               
               % If no modes exist, we delete edges.
               if isempty(relevantModes)
                    validEdgeIdx(edgeItr) = 0;
                    continue;
               end
               
               % Go through every relevant mode, and get the most likely
               % distribution's id.
               probs = zeros(size(relevantModes,1),1);
               for relevantModeItr = 1:size(relevantModes,1)
                    try
                         probs(relevantModeItr) = mvnpdf(edgeCoords(edges(edgeItr,3), :), relevantModes(relevantModeItr,4:5), [relevantModes(relevantModeItr,6:7); relevantModes(relevantModeItr,8:9)]);
                    catch %#ok<CTCH>
                         probs(relevantModeItr) = mvnpdf(edgeCoords(edges(edgeItr,3), :), relevantModes(relevantModeItr,4:5), relevantModes(relevantModeItr,[6, 9]));
                    end
               end
               
               % Weight probability densities by modes.
               probs = probs .* relevantModes(relevantModeItr,12);
               [~, newLabelIdx] = max(probs);
               newLabel = int32(relevantModes(newLabelIdx(1),3));
               edges(edgeItr,3) = newLabel;
          end
          edges = edges(validEdgeIdx, :);
          graphLevel(graphLevelItr).adjInfo = edges;
     end
end

