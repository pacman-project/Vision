%> Name: applyMove
%>
%> Description: This function applies the defined move to experts and
%> updates all data structures to reflect this change.
%>
%> @param 
%>  
%> @retval 
%>
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 12.08.2015
function [ tempExperts, tempParseTrees, tempNodeIds, tempNodeCoords, tempOrNodeChoices, tempOrNodeChoiceCounts ] = ...
     applyMove(expertToMove, parentExpertCenter, expertChildren, experts, parseTrees, nodeIds, nodeCoords, ...
          move, nodeAngles, tempOrNodeChoices, tempOrNodeChoiceCounts, nodes, vocabulary, samplingMethod, numberOfAngles, curItr)
     
     % Set step size.
     stepSize = 1;

     % Apply all moves.
     switch move(1)
          case 1
               % Assign initial stuff.
               tempExperts = experts;
               tempParseTrees = parseTrees;
               tempNodeIds = nodeIds;
               tempNodeCoords = nodeCoords;
               
               %% Position move. We'll move sub-parts around.
               % First, we obtain children positions.
              expertChildrenCoords = nodeCoords(expertChildren, :);
              for expertChildItr = 1:numel(expertChildren)
                   offsets = getMove(move(expertChildItr + 1), stepSize);
                   expertChildrenCoords(expertChildItr, :) = expertChildrenCoords(expertChildItr, :) + offsets;
              end

              % Find central point, and move children to so that their
              % center stays the same.
              centerPoint = round((max(expertChildrenCoords, [], 1) + min(expertChildrenCoords, [], 1))/2);
              expertChildrenCoords = expertChildrenCoords + repmat((parentExpertCenter - centerPoint), size(expertChildrenCoords,1),1);
              nodePosDiffs = expertChildrenCoords - nodeCoords(expertChildren, :);

              % Finally, update children (and their entire parse
              % trees) to reflect new positions.
              for expertChildItr = 1:numel(expertChildren)
                   idx = parseTrees(:,curItr+1) == expertChildren(expertChildItr);
                   tempExperts(idx, 2:3) = ...
                        tempExperts(idx, 2:3) + repmat(nodePosDiffs(expertChildItr,:), nnz(idx), 1);
                   movedExperts = unique(parseTrees(idx, curItr+1:end));
                   tempNodeCoords(movedExperts, :) = tempNodeCoords(movedExperts, :) + ...
                        repmat(nodePosDiffs(expertChildItr,:), numel(movedExperts), 1);
              end
          case 2
               %% OR node change move.
               orNodeChanges = [expertToMove, move(2)];
               [ tempExperts, tempParseTrees, tempNodeIds, tempNodeCoords, tempOrNodeChoices, tempOrNodeChoiceCounts ] = projectNode( nodes, vocabulary, samplingMethod, orNodeChanges );
          case 3
               %% Rotation move.
               
     end
end



function offsets = getMove(move, stepSize)
    switch move
          case 1
               offsets = [-stepSize, 0];
          case 2
               offsets = [stepSize, 0];
          case 3
               offsets = [0, -stepSize];
          case 4
               offsets = [0, stepSize];
%          case 5
%                offsets = [stepSize, stepSize];
%          case 6 
%                offsets = [stepSize, -stepSize];
%          case 7
%                offsets = [-stepSize, stepSize];
%          case 8
%                offsets = [-stepSize, -stepSize];
%          case 9
%                offsets = [0,0];
         case 5
               offsets = [0,0];
     end
end
