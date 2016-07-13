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
function [ graphLevel ] = assignEdgeLabels(graphLevel, modes, modeProbArr, edgeCoords, levelItr, debugFolder)
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
     orgEdgeCount = size(cat(1, allEdges{:}),1);
     
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
    secEdgeCount = size(cat(1, assignedEdges{:}),1);
     
    display([num2str(orgEdgeCount - secEdgeCount) ' (%' num2str(100 * (secEdgeCount - orgEdgeCount)/orgEdgeCount) ') edges have been removed since collected statistics are too tight.']);
    display('Balancing edges further so only doubly-linked edges remain.');
    
    %% Remove single-link edges.
    adjacentNodes = assignedEdges;
    validIdx = cellfun(@(x) ~isempty(x), adjacentNodes);
    adjacentNodes(validIdx) = cellfun(@(x) x(:,2), adjacentNodes(validIdx), 'UniformOutput', false);
    parfor graphLevelItr = 1:numel(graphLevel)
        curEdges = assignedEdges{graphLevelItr};
        curAdjacentNodes = adjacentNodes{graphLevelItr};
        if ~isempty(curAdjacentNodes)
            validEdgeIdx = ones(numel(curAdjacentNodes),1) > 0;
            for secItr = 1:numel(validEdgeIdx)
                if ~ismembc(int32(graphLevelItr), adjacentNodes{curAdjacentNodes(secItr)})
                    validEdgeIdx(secItr) = 0;
                end
            end
            if nnz(validEdgeIdx) ~= numel(validEdgeIdx)
                assignedEdges(graphLevelItr) = {curEdges(validEdgeIdx, :)};
            end
        end
    end
    finalEdgeCount = size(cat(1, assignedEdges{:}),1);
    display(['A further ' num2str(secEdgeCount - finalEdgeCount) ' (%' num2str(100 * (secEdgeCount - finalEdgeCount)/orgEdgeCount) ') edges have been removed to ensure double linking.']);
    
    edgeCounts = cellfun(@(x) size(x,1), assignedEdges);
    if usejava('jvm')
        figure('Visible', 'off'), hist(edgeCounts, 0:max(edgeCounts));
        saveas(gcf, [debugFolder '/level' num2str(levelItr) 'EdgHist_AfterRemoval.png']); 
    end
    
    % Assign final edges.
    [graphLevel.adjInfo] = deal(assignedEdges{:});
    
     % Report on deleted edges.
%     display(['%' num2str(round(100*(validEdgeCount/totalEdgeCount)))...
%         ' of all edges have been preserved. Total: ' num2str(totalEdgeCount) ...
%         ' edges. Deleted:' num2str(totalEdgeCount - validEdgeCount) ' edges.']); 
     
     clearvars -except graphLevel
end