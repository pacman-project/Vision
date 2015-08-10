function [ graphLevel ] = assignEdgeLabels( graphLevel, modes, edgeCoords)
     allEdges = cat(1, graphLevel.adjInfo);
     nodeIds = [graphLevel.labelId]';
     
     if isempty(allEdges)
          return;
     end
     
     % Get node ids of edges.
     for graphLevelItr = 1:numel(graphLevel)
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
               [~, newLabelIdx] = max(probs);
               newLabel = int32(relevantModes(newLabelIdx(1),3));
               edges(edgeItr,3) = newLabel;
          end
          edges = edges(validEdgeIdx, :);
          graphLevel(graphLevelItr).adjInfo = edges;
     end
end

