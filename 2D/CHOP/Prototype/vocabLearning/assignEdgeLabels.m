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
     
     % Get node ids of edges.
     parfor graphLevelItr = 1:numel(graphLevel)
          edges = graphLevel(graphLevelItr).adjInfo;

          % Edges empty, do nothing.
          if isempty(edges)
               continue;
          end
          
          edgeNodeLabels = nodeIds(edges(:, 1:2)); %#ok<PFBNS>
          if size(edgeNodeLabels,2) ~= 2
               edgeNodeLabels = edgeNodeLabels';
          end
          
          % Assign each edge a to a relevant mode.
          probArr = zeros(size(edges,1),1, 'single');
          graphLevel(graphLevelItr).realEdgeLabels = edges(:, 3);
          for edgeItr = 1:size(edges,1)
               relevantIdx = modes(:, 1) == edgeNodeLabels(edgeItr,1) & modes(:,2) == edgeNodeLabels(edgeItr,2); %#ok<PFBNS>
               relevantModes = modes(relevantIdx,:);
               relevantCoords = edgeCoords(edges(edgeItr,3),:); %#ok<PFBNS>
               clusterProbs = modeProbArr(relevantIdx,relevantCoords(1),relevantCoords(2)); %#ok<PFBNS>
               [probability, clusterId] = max(clusterProbs);
               probArr(edgeItr) = probability;
               newLabel = int32(relevantModes(clusterId,3));
               edges(edgeItr,3) = newLabel;
          end
          graphLevel(graphLevelItr).adjInfo = edges;
          graphLevel(graphLevelItr).edgeProbabilities = probArr;
     end
end