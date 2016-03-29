%% Function to generate moves for optimization process.
%> Name: generateMoves
%>
%> Description: This function generates a set of moves for SGD-based
%> optimization.
%>
%> @param stochasticMoves Number of moves
%> @param numberOfChildren Number of children to be moved around.
%> @moveFlags Allowed moves (3x1 array). moveFlags(1) refers to positional
%> moves. moveFlags(2) is for OR node moves. moveFlags(3) is for rotational
%> moves.
%> @param expertOrNodeChoice The OR node choice that's already made.
%> @param expertOrNodeChoiceCount Total number of OR node choices.
%> 
%> @retval moves Randomly generated unique moves.
%>  
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 24.03.2016
function moves = generateMoves(stochasticMoves, numberOfChildren, moveFlags, expertOrNodeChoice, expertOrNodeChoiceCount, topLevelFlag)
     generatedMoveCount = stochasticMoves * 2;
     moves = zeros(generatedMoveCount, numberOfChildren + 1);
     
     % If there's only 1 type of OR node choice, we set the relevant flag
     % to 0. If only this mode is active, we return.
     if expertOrNodeChoiceCount == 1 && moveFlags(2) == 1 && nnz(moveFlags) == 1
          if nnz(moveFlags) == 1
               moves = [];
               return;
          else
               moveFlags(2) = 0;
          end
     end
     
     % If this is a top node with nothing around, we do not perform OR node
     % or rotation operations.
     if topLevelFlag
         moveFlags(2:3) = 0;
         if ~moveFlags(1)
             moves = [];
         end
     end
     
     % Finally, get all available types.
     availableTypes = find(moveFlags);
     
    % If there's only one type of move allowed, we don't need to perform
    % initial selection.
    if nnz(moveFlags) == 1
         moves(:,1) = availableTypes;
    else
         moves(:,1) = ceil(rand(generatedMoveCount, 1) * nnz(moveFlags));
         moves(moves(:,1)<1,1) = 1;
         moves(:,1) = availableTypes(moves(:,1));
    end

    % 1) Now, let's also generate position moves for each entry, where needed.
    numberOfPosMoves = nnz(moves(:,1) == 1);
    if numberOfPosMoves > 0
        posMoves = ceil(rand(numberOfPosMoves, numberOfChildren) * 5);
        posMoves(posMoves < 1) = 1;
        moves(moves(:,1) == 1, 2:end) = posMoves;
    end
    
    % 2) Finally, we implement OR nodes.
    numberOfOrMoves = nnz(moves(:,1) == 2);
    orNodeOptions = setdiff(1:expertOrNodeChoiceCount, expertOrNodeChoice);
    if numberOfOrMoves > 0
        orMoves = ceil(rand(numberOfOrMoves, 1) * (expertOrNodeChoiceCount-1));
        orMoves(orMoves < 1) = 1;
        moves(moves(:,1) == 2, 2) = orNodeOptions(orMoves);
    end
    
    % 3) Rotation entries are needed.
    numberOfRotMoves = nnz(moves(:,1) == 3);
    if numberOfRotMoves > 0
        rotMoves = ceil(rand(numberOfRotMoves, 1) * 2);
        rotMoves(rotMoves < 1) = 1;
        moves(moves(:,1) == 3, 2) = rotMoves;
    end
    
   % Remove invalid position moves.
   if topLevelFlag
      moves = moves(~(moves(:,1) == 1 & range(moves(:,2:end),2) == 0), :);  
   end
   
    % The moves are now generated. Here, we remove duplicate rows.
    moves = unique(moves, 'stable', 'rows');
    
    % If we have too many moves, trim the matrix.
    if size(moves,1) > stochasticMoves
         moves = moves(1:stochasticMoves,:);
    end
end