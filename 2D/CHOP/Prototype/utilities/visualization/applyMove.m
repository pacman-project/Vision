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
function [ newExperts, newSubChildrenExperts, expertChildren, newAngle ]  = ...
                        applyMove(nodeCoords, newExperts, expertChildren, newSubChildrenExperts, move, nodeAngle, numberOfRealFilters, numberOfFilters)
     
     % Set step size.
     stepSize = 1;
     newAngle = nodeAngle;

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
                    filterIds = round(((180/numberOfRealFilters) * (0:(numberOfRealFilters-1))) / (180/numberOfFilters))' + 1;
                    newExperts(:,1) = filterIds(newExperts(:,1));
              end
          case 2
               %% Or node change. We don't need to do anything here.
               
          case 3
               %% Rotation move.
               if move(2) == 1
                    newAngle = mod(nodeAngle - 1, 2 * numberOfFilters);
               else
                    newAngle = mod(nodeAngle + 1, 2 * numberOfFilters);
               end
     end
     
     % Rendering part! We have performed the required operations, now let's
     % obtain the experts and rotate everything around the center point.
     if newAngle > 0
          % Create data structures.
          newExperts = double(newExperts);
          expertAngles = mod(newAngle + (newExperts(:,1) - 1), numberOfFilters) + 1;
          numberOfExperts = size(newExperts,1);
          expertCoords = newExperts(:,2:3);
          
          % Rotation matrices
          center = repmat(nodeCoords, numberOfExperts, 1);
          theta = -(2 * pi * newAngle / (numberOfFilters));
          
          % Perform rotation.
          R = [cos(theta) -sin(theta); sin(theta) cos(theta)]; 
          expertCoords = (R * ((expertCoords - center))' + center')';
          
          % Get new experts.
          newExperts = round([expertAngles, expertCoords]);
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
