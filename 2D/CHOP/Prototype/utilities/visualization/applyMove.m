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
function [ newExperts, newSubChildrenExperts, expertChildren, rotatedExpertChildren, newOrNodeChoice, newAngle ]  = ...
                        applyMove(nodeCoords, newExperts, expertChildren, newSubChildrenExperts, move, expertOrNodeChoice, nodeAngle, numberOfRealFilters, numberOfFilters)
     
     % Set step size.
     stepSize = unidrnd(3);
     newAngle = nodeAngle;
     newOrNodeChoice = expertOrNodeChoice;
     nodeCoords = double(nodeCoords);
     filterIds = round(((180/numberOfRealFilters) * (0:(numberOfRealFilters-1))) / (180/numberOfFilters))' + 1;

     % Apply all moves.
     switch move(1)
          case 1
               %% Position move. We'll move sub-parts around.
               % First, we obtain children positions.
              expertChildrenCoords = double(expertChildren(:,2:3));
              for expertChildItr = 1:size(expertChildren,1)
                   offsets = getMove(move(expertChildItr + 1), stepSize);
                   expertChildrenCoords(expertChildItr, :) = expertChildrenCoords(expertChildItr, :) + offsets;
              end

              % Find central point, and move children to so that their
              % center stays the same.
              centerPoint = round((max(expertChildrenCoords, [], 1) + min(expertChildrenCoords, [], 1))/2);
              expertChildrenCoords = expertChildrenCoords + repmat((nodeCoords - centerPoint), size(expertChildrenCoords,1),1);
              nodePosDiffs = expertChildrenCoords - double(expertChildren(:,2:3));
              
              % Update the children.
              expertChildren(:,2:3) = int32(expertChildrenCoords);
              for nodeItr = 1:size(expertChildren,1)
                   tempVar = double(newSubChildrenExperts{nodeItr});
                   tempVar(:,2:3) = tempVar(:,2:3) + repmat(nodePosDiffs(nodeItr,:), size(tempVar,1), 1);
                   newSubChildrenExperts{nodeItr} = int32(tempVar);
              end
              newExperts = cat(1, newSubChildrenExperts{:});
              
              % Rotation ids
              if numberOfRealFilters < numberOfFilters
                    newExperts(:,1) = filterIds(newExperts(:,1));
              end
          case 2
               %% Or node change. We don't need to do anything here.
               newOrNodeChoice = move(2);
          case 3
               %% Rotation move.
               if move(2) == 1
                    newAngle = mod(nodeAngle - 1, 2 * numberOfFilters);
               else
                    newAngle = mod(nodeAngle + 1, 2 * numberOfFilters);
               end
          case 4
               %% Only for rendering.
     end
     
     % Save rotated expert children.
     rotatedExpertChildren = double(expertChildren);
     
     % Rendering part! We have performed the required operations, now let's
     % obtain the experts and rotate everything around the center point.
     if newAngle ~= 0 
          % Create data structures.
          expertChildren = double(expertChildren);
          center = repmat(nodeCoords, size(expertChildren, 1), 1);
          theta = -(pi * newAngle / (numberOfFilters));
          R = [cos(theta) -sin(theta); sin(theta) cos(theta)]; 
          
          % Rotate children.
          rotatedExpertChildren(:,2:3) = round((R * ((expertChildren(:,2:3) - center))' + center')');
          childrenOffsets = int32(rotatedExpertChildren - expertChildren);
          
          % Rotate low-level experts of children. First, we shift their
          % locations to the new places.
          newExperts = newSubChildrenExperts;
          for childItr = 1:size(rotatedExpertChildren,1)
              tempExperts = newSubChildrenExperts{childItr};
              tempExperts(:,2:3) = tempExperts(:,2:3) + repmat(childrenOffsets(childItr,2:3), size(tempExperts,1), 1);
              newSubChildrenExperts{childItr} = tempExperts;
              
              % Rotate low-level experts.
              center = repmat(rotatedExpertChildren(childItr,2:3), size(tempExperts,1), 1);
              tempExperts(:,2:3) = round((R * ((double(tempExperts(:,2:3)) - center))' + center')');
              
              % First, assign crude ids.
              if numberOfRealFilters < numberOfFilters
                    tempExperts(:,1) = filterIds(tempExperts(:,1));
              end
              
              % Then, update the filter ids (for finer detail)
              tempExperts(:,1) = mod(newAngle + (tempExperts(:,1) - 1), numberOfFilters) + 1;
              
              % Save final (rendered) experts.
              newExperts{childItr} = tempExperts;
          end
          newExperts = cat(1, newExperts{:});
          expertChildren = int32(expertChildren);
          rotatedExpertChildren = int32(rotatedExpertChildren);
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
         case 5
               offsets = [0,0];
     end
end
