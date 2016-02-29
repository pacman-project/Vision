%> Name: optimizePoE
%>
%> Description: We project high level nodes to an image, calculate
%> gradients for nodes across all layers, and slowly move them around to
%> climb the data likelihood function. 
%>
%> @param nodes High-level nodes in format of (realLabelId, posX, posY,
%levelItr; ...]. 
%> @param vocabulary Part vocabulary.
%> @param imageSize Image size in pixels.
%> @param options Program options.
%> 
%> @retval optimizedImg Optimized image (backprojected).
%>
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 07.02.2016
function [ modalImg, experts ] = optimizeImagination( nodes, vocabulary, imageSize, prevImageSize, filters, visFilters, sampleItr, batchFlag, datasetName, fileName)
     stepSize = 1;
     minPixelValue = 1/255;
     stopVal = 5;
     if ~batchFlag
          maxSteps = 5 * size(nodes,1);
     else
          maxSteps = 5;
     end
     
     % If the file name is given, use that one.
     if nargin < 9
          fileName = num2str(nodes(1,1));
     end
     folderName = [pwd '/optimization/' datasetName '/level' num2str(nodes(1,4)) '/' fileName '_sample' num2str(sampleItr)];
     
     % Create ourput folder.
     if ~exist(folderName, 'dir')
          mkdir(folderName);
     end
     
     movesPerChild = 2;
     poeFlag = true;
     positionFlag = false;
     imageChoices = prevImageSize * prevImageSize;
     for filtItr = 1:numel(visFilters)
          tempFilt = visFilters{filtItr};
          tempFilt = tempFilt/max(max(tempFilt));
          tempFilt(tempFilt<minPixelValue) = minPixelValue;
          visFilters{filtItr} = tempFilt;
     end
%     filters = visFilters;

     % Shut down warnings.
     warning('off','all');
     
     %% First, imagine layer 1 nodes.
     [experts, parseTrees, nodeIds, nodeCoords] = projectNode(nodes, vocabulary, 1, 'modal');
     experts = double(experts);
     
     % To work with more precise coordinates, we use imprecise coordinates
     % here.
     experts(:,2:3) = nodeCoords(unique(parseTrees(:,end)),:);
     
    % Center the nodes.
    minX = min(experts(:,2));
    maxX = max(experts(:,2));
    minY = min(experts(:,3));
    maxY = max(experts(:,3));
    midPoint = round([(minX+maxX)/2, (minY+maxY)/2]);
    experts(:,2:3) = (experts(:,2:3) - repmat(midPoint, size(experts,1), 1))+ round(repmat(imageSize/2, size(experts,1), 1));
    nodeCoords = (nodeCoords - repmat(midPoint, size(nodeCoords,1), 1) + round(repmat(imageSize/2, size(nodeCoords,1), 1)));
    
    %% Second, we start a gradient-descent procedure to move things around so they match.
    modalImg = [];
    varMat = [];
    likelihoodMat = [];
     
    % Only move nodes from level below.
    steps = 0;
    curItr = 1;
    curLevelItr = (size(parseTrees,2) - curItr) + 1;
    currentLevel = vocabulary{curLevelItr};
    imageArr = cell(maxSteps,1);
    posLikelihoodArr = zeros(maxSteps,1);
    poeLikelihoodArr = zeros(maxSteps,1);
    diffImageArr = cell(maxSteps,1);
    diffImageArr{1} = zeros(imageSize);
    
    % Save 
    batchExperts = experts;
    batchNodeCoords = nodeCoords;
    
    %% Continue with gradient descent until optimized.
    while steps < maxSteps && curItr <= size(parseTrees,2)-1
         
         % Select experts to move.
         [modifiedExperts, ~] = unique(parseTrees(:,curItr), 'first');
         modifiedExperts = modifiedExperts';
         
         %% Batch gradient descent requires likelihood calculations here.
         % Get reference likelihood value.       
         if batchFlag
              [refModalImg, ~, ~, ~] = obtainPoE(experts, modalImg, varMat, likelihoodMat, imageSize, visFilters);
              if ~poeFlag
                   refLikelihoodVal = 0;
              else
                   [~, ~, refLikelihoodMat, refLikelihoodVal] = obtainPoE(experts, modalImg, varMat, likelihoodMat, imageSize, filters);
                   % Normalize by number of choices.
       %            refLikelihoodVal = refLikelihoodVal / imageChoices;
              end
              refModalMask = refModalImg > minPixelValue;

              % Print images.
              imgName = [folderName '/' num2str(nodes(1,1)) '_sample' num2str(sampleItr) '_step' num2str(steps) '_level' num2str(curLevelItr) '.png'];
              maskName = [folderName '/mask_' num2str(nodes(1,1)) '_sample' num2str(sampleItr) '_step' num2str(steps) '_level' num2str(curLevelItr) '_Mask.png'];
              imwrite(refModalImg / max(max(refModalImg)), imgName);
              imwrite(refModalMask, maskName);

              % Save image.
              imageArr{steps+1} = refModalImg / max(max(refModalImg));
              poeLikelihoodArr(steps+1) = refLikelihoodVal;
         end
         
         
         %% Here, we create random perturbations in the joint distribution
         % space for each expert. Then, we'll evaluate the coupled gradient 
         % ( position gradients + Product of Experts gradients). Then,
         % we'll make the best move based on this joint gradient.
         % Just in case, we sample two times the normal amount not to have
         % repetitions.  
         moveFlag = false;
         for expertItr = 1:numel(modifiedExperts)
              % Get reference likelihood value.       
              if ~batchFlag
                   [refModalImg, ~, ~, ~] = obtainPoE(experts, modalImg, varMat, likelihoodMat, imageSize, visFilters);
                   if ~poeFlag
                        refLikelihoodVal = 0;
                   else
                        [~, ~, refLikelihoodMat, refLikelihoodVal] = obtainPoE(experts, modalImg, varMat, likelihoodMat, imageSize, filters);
                        % Normalize by number of choices.
            %            refLikelihoodVal = refLikelihoodVal / imageChoices;
                   end
                   refModalMask = refModalImg > minPixelValue;

                   % Print images.
                   imgName = [folderName '/' num2str(nodes(1,1)) '_sample' num2str(sampleItr) '_step' num2str(steps) '_level' num2str(curLevelItr) '.png'];
                   maskName = [folderName '/mask_' num2str(nodes(1,1)) '_sample' num2str(sampleItr) '_step' num2str(steps) '_level' num2str(curLevelItr) '_Mask.png'];
                   imwrite(refModalImg / max(max(refModalImg)), imgName);
                   imwrite(refModalMask, maskName);

                   % Save image.
                   imageArr{steps+1} = refModalImg / max(max(refModalImg));
                   poeLikelihoodArr(steps+1) = refLikelihoodVal;
              end
              
              % Mark expert as immobile. Then, we'll see if we can move it
              % to obtain a lower likelihood.
              expertMoveFlag = false;
              expertLabelId = nodeIds(modifiedExperts(expertItr));
              expertNode = currentLevel(expertLabelId);
              
               % First of all, we select the relevant combination of nodes.
               expertChildren = unique(parseTrees(parseTrees(:, curItr) == modifiedExperts(expertItr), curItr+1), 'first');
               childrenLabelCombination = nodeIds(expertChildren);
               relevantCombination = find(ismember(childrenLabelCombination, expertNode.childrenLabelDistributions(:,1:(end-1)), 'rows'));
               
               % Calculate pos probability denominator.
               posProbDenom  = imageChoices ^ (numel(expertChildren));

               % Then, we obtain the position pdfs.
               childrenPosDistributions = expertNode.childrenPosDistributions{1}{relevantCombination}; %#ok<FNDSB>
              
              % Next thing we do, we sample from the available moves
              % (joint).
              stochasticMoves = round(movesPerChild^numel(expertChildren));
              moves = ceil(rand(stochasticMoves*5, numel(expertNode.children)) * 5);
              moves(moves<1) = 1;
              moves = unique(moves, 'stable', 'rows');
              if numel(modifiedExperts) == 1 && curItr == 1
                   moves = moves(range(moves,2) ~= 0, :);
              else
                   moves = moves(~(range(moves,2) == 0 & moves(:,1) == 5), :);
              end
              if size(moves,1) > stochasticMoves
                   moves = moves(1:stochasticMoves,:);
              end
              numberOfMoves = size(moves,1);
              
              %% Calculate existing prediction's position likelihood.
              childExpertCenter = expertChildren(1);
              if numel(expertChildren) > 1
                   childExpertPers = expertChildren(2:end);
              else
                   childExpertPers = [];
              end

              % If positions are involved, 
              if positionFlag
                   expertPosLikelihoodVal = GetTreeLikelihood(childExpertCenter, childExpertPers, nodeCoords, childrenPosDistributions, posProbDenom);
              else
                   expertPosLikelihoodVal = 0;
              end
              
              %% Update expert positions based on the perturbations.
              maxExperts = experts;
              maxNodeCoords = nodeCoords;
              gradients = zeros(numberOfMoves,1);
              parentExpertCenter = nodeCoords(modifiedExperts(expertItr), :);
              maxMove = 0;
              maxNewPoELikelihoodVal = 0;
              maxLikelihoodMat = [];
              maxNewPosLikelihoodVal = 0;
              
              % Before we start moving, we calculate an initial product of
              % experts by not taking the expert in question into account.
              if poeFlag              
                   oldExperts = experts(~(parseTrees(:, curItr) == modifiedExperts(expertItr)),:);
                   if ~isempty(oldExperts)
                        [preModalImg, preVarMat, preLikelihoodMat, ~] = obtainPoE(oldExperts, modalImg, varMat, likelihoodMat, imageSize, filters);
                   else
                        preModalImg = [];
                        preVarMat = [];
                        preLikelihoodMat = [];
                   end
              end
              
              addedExperts = parseTrees(:, curItr) == modifiedExperts(expertItr);
              
              % Calculate moves, and their gradients.
              for moveItr = 1:numberOfMoves
                   %% Perturb the experts based on assigned moves.
                   tempExperts = experts;
                   
                   % First, we obtain children positions.
                   expertChildrenCoords = nodeCoords(expertChildren, :);
                   tempNodeCoords = nodeCoords;
                   for expertChildItr = 1:numel(expertChildren)
                        offsets = getMove(moves(moveItr, expertChildItr), stepSize);
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
                   
                   %% Calculate new position likelihood.
                   if positionFlag
                        newPosLikelihoodVal = GetTreeLikelihood(childExpertCenter, childExpertPers, tempNodeCoords, childrenPosDistributions, posProbDenom);
                   else
                        newPosLikelihoodVal = 0;
                   end
                   
                   %% Calculate pixel-level likelihoods
                   expertJointLikelihoodVal = refLikelihoodVal + expertPosLikelihoodVal;
                   if poeFlag
                        if ~isinf(newPosLikelihoodVal)
                             [~, ~, newLikelihoodMat, ~] = obtainPoE(tempExperts, preModalImg, preVarMat, preLikelihoodMat, imageSize, filters, addedExperts);
                             newPoELikelihoodVal = sum(sum(newLikelihoodMat));
                             newVal = newPosLikelihoodVal + newPoELikelihoodVal;
                        else
                             newVal = -Inf;
                             newPoELikelihoodVal = 0;
                        end
                   else
                        newVal = newPosLikelihoodVal;
                        newPoELikelihoodVal = 0;
                   end
                   
                   %% Finally, calculate the difference between the likelihoods and update gradients.
                   newDiffVal = newVal - expertJointLikelihoodVal;
                   if newDiffVal > stopVal && newDiffVal > max(gradients)
                        maxExperts = tempExperts;
                        maxNewPoELikelihoodVal = newPoELikelihoodVal;
                        maxNewPosLikelihoodVal = newPosLikelihoodVal;
                        maxLikelihoodMat = newLikelihoodMat;
                        maxNodeCoords = tempNodeCoords;
                        expertMoveFlag = true;
                        maxMove = moveItr;
                   end
                   
                   % Calculate gradients.
                    gradients(moveItr) = newDiffVal;
              end
              
              % If we haven't moved anything, nothing to do here.
              if ~expertMoveFlag
                   continue;
              else
                   moveFlag = true;
                   % If we're not performing batch gradient descent, let's
                   % break and update our reference likelihood values. 
                   changedExperts = unique(parseTrees(parseTrees(:,curItr) == modifiedExperts(expertItr),(curItr+1):end));
                   display(['Expert ' num2str(modifiedExperts(expertItr)) ' has been modified.']);
                   display(['Position likelihood changed from ' num2str(round(expertPosLikelihoodVal)) ' to ' num2str(round(maxNewPosLikelihoodVal)) '.']);
                   if maxNewPosLikelihoodVal > expertPosLikelihoodVal && maxNewPoELikelihoodVal > refLikelihoodVal
                        display('Both likelihoods improved!');
                   elseif maxNewPosLikelihoodVal > expertPosLikelihoodVal
                        display('Position distributions pulled the strings.');
                   elseif maxNewPoELikelihoodVal > refLikelihoodVal
                        display('Product of experts predictions dominated.');
                   end
                   if ~batchFlag
                        % This expert has been modified. Let's update the
                        % existing data structures to reflect this change.
                        display(['PoE likelihood changed from ' num2str(round(refLikelihoodVal)) ' to ' num2str(round(maxNewPoELikelihoodVal)) '.']);
                        posLikelihoodArr(steps+1) = newPosLikelihoodVal - expertPosLikelihoodVal;
                        refLikelihoodVal = maxNewPoELikelihoodVal - gradients(maxMove);
                        experts(parseTrees(:,curItr) == modifiedExperts(expertItr),:) = maxExperts(parseTrees(:,curItr) == modifiedExperts(expertItr),:);
                        nodeCoords(changedExperts,:) = maxNodeCoords(changedExperts,:);
                        steps = steps+1;
                   else
                        display(['--For this expert, PoE likelihood changed from ' num2str(round(refLikelihoodVal)) ' to ' num2str(round(maxNewPoELikelihoodVal)) '.']);
                        batchExperts(parseTrees(:,curItr) == modifiedExperts(expertItr),:) = maxExperts(parseTrees(:,curItr) == modifiedExperts(expertItr),:);
                        batchNodeCoords(changedExperts,:) = maxNodeCoords(changedExperts,:);
                   end
              end
         end
         
         % If none of the experts was able to move, we move on to the next
         % level's realizations.
         if ~moveFlag
              curItr = curItr + 1;
              curLevelItr = (size(parseTrees,2) - curItr) + 1;
              currentLevel = vocabulary{curLevelItr};
         else
              % If batch flag is active, we should calculate a final
              % likelihood. Otherwise, we can use existing values that come
              % from the halted loop. 
              if batchFlag
                   experts = batchExperts;
                   nodeCoords = batchNodeCoords;
                   [~, ~, maxLikelihoodMat, maxNewPoELikelihoodVal] = obtainPoE(experts, modalImg, varMat, likelihoodMat, imageSize, filters);
                   display(['PoE likelihood changed from ' num2str(round(refLikelihoodVal)) ' to ' num2str(round(maxNewPoELikelihoodVal)) '.']);
                   % Increase steps.
                   steps = steps + 1;              
                   
                   % Obtain the difference of the likelihood image.
                   if ~isempty(maxLikelihoodMat)     
                        diffLikelihoodMat = abs(maxLikelihoodMat - refLikelihoodMat);
                        imwrite(diffLikelihoodMat / max(max(diffLikelihoodMat)), [folderName '/diff_' num2str(nodes(1,1)) '_sample' num2str(sampleItr) '_step' num2str(steps) '_level' num2str(curLevelItr) '.png']);
                        diffImageArr{steps+2} = diffLikelihoodMat;
                   end
              end


         end
    end
    
    % Create a gif image out of diagnostic information.
    if steps>0
         if steps<maxSteps
              posLikelihoodArr = posLikelihoodArr(1:steps+1);
              poeLikelihoodArr= poeLikelihoodArr(1:steps+1);
              diffImageArr = diffImageArr(1:steps+1);
         end
         diffMax = max(cellfun(@(x) max(max(x)), diffImageArr));
         if diffMax>0
              diffImageArr = cellfun(@(x) x/diffMax, diffImageArr, 'UniformOutput', false);
         end
         posPadding = (max(posLikelihoodArr) - min(posLikelihoodArr))/4;
         posLimits = [min(posLikelihoodArr) - posPadding, max(posLikelihoodArr) + posPadding];
         if numel(unique(posLimits)) == 1
              posLimits = [posLimits(1)-1, posLimits(1)+1];
         end
         poePadding = (max(poeLikelihoodArr) - min(poeLikelihoodArr))/4;
         poeLimits = [min(poeLikelihoodArr) - poePadding, max(poeLikelihoodArr) + poePadding];
         if numel(unique(poeLimits)) == 1
              poeLimits = [poeLimits(1)-1, poeLimits(1) +1];
         end
         fileName = [folderName '.gif'];
         for stepItr = 1:min(maxSteps, (steps+1))
              figure('Visible', 'off'), hold on;
              axis square
              subplot(2,2,1), imshow(imageArr{stepItr});
              title('Imagined Data')

              subplot(2,2,2), imshow(diffImageArr{stepItr});
              title('Likelihood differences')
              
              subplot(2,2,3), plot(1:stepItr, posLikelihoodArr(1:stepItr));
              ylim(posLimits);
              if stepItr>1
                    xlim([1, min(maxSteps, (steps+1))]);
              end
              title('Change in position likelihood')
              
              subplot(2,2,4), plot(1:stepItr, poeLikelihoodArr(1:stepItr));
              ylim(poeLimits);
              if stepItr>1
                    xlim([1, min(maxSteps, (steps+1))]);
              end
              title('Product of experts likelihood')
              
              hold off;
              saveas(gcf,  [folderName '_temp.png']);
              im=imread([folderName '_temp.png']);
              [imind,cm] = rgb2ind(im, 256);
              
              if stepItr == 1
                   imwrite(imind, cm, fileName, 'LoopCount', inf, 'DelayTime',2);
              else
                   imwrite(imind, cm, fileName, 'WriteMode', 'append');
              end
              close(gcf);
         end
    end
    
    % Turn warnings on.
     warning('on','all');
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


%% Function to calculate parse tree likelihood. 
function [expertPosLikelihoodVal] = GetTreeLikelihood(childExpertCenter, childExpertPers, nodeCoords, childrenPosDistributions, posProbDenom)
    if ~isempty(childExpertPers)
         sampledPos = nodeCoords(childExpertPers,:) - repmat(nodeCoords(childExpertCenter,:), numel(childExpertPers),1);
         refSampledPos = [];
         for childItr = 1:size(sampledPos,1)
              refSampledPos = [refSampledPos, sampledPos(childItr,:)]; %#ok<AGROW>
         end
         try
              posProbability = pdf(childrenPosDistributions, refSampledPos)  / posProbDenom;
         catch %#ok<CTCH>
              posProbability = 0;
         end
         expertPosLikelihoodVal = log(posProbability);
    else
         expertPosLikelihoodVal = 0;
    end
end
