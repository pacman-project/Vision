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
function [ modalImg, experts ] = optimizeImagination( nodes, vocabulary, imageSize, prevImageSize, filters, visFilters, sampleItr, batchFlag, datasetName, likelihoodLookupTable, fileName)
     stepSize = 1;
     stopVal = 1;
     if ~batchFlag
          maxSteps = 5 * size(nodes,1);
     else
          maxSteps = 5;
     end
     
     % If the file name is given, use that one.
     sampleString = 'modal';
     if nargin < 11
          fileName = num2str(nodes(1,1));
          sampleString = 'sample';
     end
     folderName = [pwd '/optimization/' datasetName '/level' num2str(nodes(1,4)) '/' fileName '_sample' num2str(sampleItr)];
     
     % Create a circular structure to move around the image.
     filterSize = size(filters,1);
     maxInfluenceRadius = ceil(filterSize/2);
     cx=maxInfluenceRadius+1;cy=cx;ix=2*maxInfluenceRadius+1;iy=ix;r=maxInfluenceRadius;
     [x,y]=meshgrid(-(cx-1):(ix-cx),-(cy-1):(iy-cy));
     c_mask=((x.^2+y.^2)<r^2);
     
     % Create ourput folder.
     if ~exist(folderName, 'dir')
          mkdir(folderName);
     end
     movesPerChild = 2;
     maxMoves = 20;
     poeFlag = true;
     positionFlag = false;
     imageChoices = prevImageSize * prevImageSize;
     
     % Shut down warnings.
     warning('off','all');
     
     %% First, imagine layer 1 nodes.
     [experts, parseTrees, nodeIds, nodeCoords] = projectNode(nodes, vocabulary, 1, sampleString);
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
    refLikelihoodMat = zeros(imageSize);
    
    % Save 
    batchExperts = experts;
    batchNodeCoords = nodeCoords;
    
    %% Continue with gradient descent until optimized.
    while steps < maxSteps && curItr <= size(parseTrees,2)-1
         
         % Select experts to move.
         [modifiedExperts, ~] = unique(parseTrees(:,curItr), 'first');
         modifiedExperts = modifiedExperts';
         
         %% We obtain the reference modal image and likelihood mat, in order to be able to select parts based on their local likelihood values.
         [refModalImg, ~, ~] = obtainPoE(experts, modalImg, likelihoodMat, imageSize, visFilters, []);
         [~, refLikelihoodMat, refLikelihoodVal] = obtainPoE(experts, modalImg, likelihoodMat, imageSize, filters, likelihoodLookupTable);
         
         % Check which experts are to be moved. We'll move the expert which
         % has lowest average likelihood.
         avgLikelihoods = zeros(numel(modifiedExperts),1);
         expertIdx = cell(numel(modifiedExperts),1);
         tempMatrix = zeros(imageSize) > 0;
         ind = sub2ind(imageSize, experts(:,2), experts(:,3));
         for expertItr = 1:numel(modifiedExperts)
              lowLevelExperts = ind(parseTrees(:,curItr) == modifiedExperts(expertItr));
              tempMatrix(lowLevelExperts) = 1;
              tempOverlapMatrix = imdilate(tempMatrix, c_mask);
              expertIdx(expertItr) = {find(tempMatrix)};
              avgLikelihoods(expertItr) = mean(refLikelihoodMat(tempOverlapMatrix));
         end
         moveFlagArr = ones(numel(modifiedExperts),1) > 0;
         
         while steps<maxSteps
              % Based on the list of matrices at hand, let's iterate over the
              % list and start fixing things by moving them around.
              [sortedExpertLikelihoods, sortedExpertIdx] = sort(avgLikelihoods);
              
              % If no experts could be moved, move on.
              if nnz(moveFlagArr) == 0
                   break;
              end
              
              % Select the expert with minimum average likelihood.
              curExpertItr = find(moveFlagArr, 1, 'first');
              tempIdx = sortedExpertIdx(find(moveFlagArr, 1, 'first'));
              
              % Move the expert!
              expertToMove = modifiedExperts(tempIdx);
              expertLabelId = nodeIds(expertToMove);
              expertNode = currentLevel(expertLabelId);
              
               % First of all, we select the relevant combination of nodes.
               expertChildren = unique(parseTrees(parseTrees(:, curItr) == expertToMove, curItr+1), 'first');
               childrenLabelCombination = nodeIds(expertChildren);
               relevantCombination = find(ismember(childrenLabelCombination, expertNode.childrenLabelDistributions(:,1:(end-1)), 'rows'));
               
               % Calculate pos probability denominator.
               posProbDenom  = imageChoices ^ (numel(expertChildren));

               % Then, we obtain the position pdfs.
               childrenPosDistributions = expertNode.childrenPosDistributions{1}{relevantCombination}; %#ok<FNDSB>
               
              % Next thing we do, we sample from the available moves
              % (joint).
              stochasticMoves = round(min(movesPerChild^numel(expertChildren), maxMoves));
              moves = ceil(rand(stochasticMoves, numel(expertNode.children)) * 5);
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
              gradients = zeros(numberOfMoves,1);
              parentExpertCenter = nodeCoords(expertToMove, :);
              
              % Before we start moving, we calculate an initial product of
              % experts by not taking the expert in question into account.
              if poeFlag              
                   oldExperts = experts(~(parseTrees(:, curItr) == expertToMove),:);
                   if ~isempty(oldExperts)
                        [preModalImg, preLikelihoodMat, ~] = obtainPoE(oldExperts, modalImg, likelihoodMat, imageSize, filters, likelihoodLookupTable);
                   else
                        preModalImg = [];
                        preLikelihoodMat = [];
                   end
              end
              addedExperts = double(parseTrees(:, curItr) == expertToMove);
              
              % Allocate space for max experts.
              maxExperts = cell(numberOfMoves,1);
              maxNewPoELikelihoodVal = zeros(numberOfMoves,1);
              maxNewPosLikelihoodVal = zeros(numberOfMoves,1);
              maxNodeCoords = cell(numberOfMoves,1);
              
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
                             [~, newLikelihoodMat, ~] = obtainPoE(tempExperts, preModalImg, preLikelihoodMat, imageSize, filters, likelihoodLookupTable, addedExperts);
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

                   % Save the info.
                   if newDiffVal > stopVal
                        maxExperts(moveItr) = {tempExperts};
                        maxNewPoELikelihoodVal(moveItr) = newPoELikelihoodVal;
                        maxNewPosLikelihoodVal(moveItr) = newPosLikelihoodVal;
                        maxNodeCoords(moveItr) = {tempNodeCoords};
                   end
                   
                   % Calculate gradients.
                    gradients(moveItr) = newDiffVal;
              end
              
              if nnz(gradients) == 0
                   moveFlagArr(curExpertItr) = 0;
                   continue;
              end
              
              % We have moved the expert, now let's update the data
              % structures.
              [~, maxMove] = max(gradients);
              
              % If we're not performing batch gradient descent, let's
              % break and update our reference likelihood values. 
              changedExperts = unique(parseTrees(parseTrees(:,curItr) == expertToMove,(curItr+1):end));
              display(['Expert ' num2str(expertToMove) ' has been modified.']);
              display(['Position likelihood changed from ' num2str(round(expertPosLikelihoodVal)) ' to ' num2str(round(maxNewPosLikelihoodVal(maxMove))) '.']);
              if maxNewPosLikelihoodVal(maxMove) > expertPosLikelihoodVal && maxNewPoELikelihoodVal(maxMove) > refLikelihoodVal
                   display('Both likelihoods improved!');
              elseif maxNewPosLikelihoodVal(maxMove) > expertPosLikelihoodVal
                   display('Position distributions pulled the strings.');
              elseif maxNewPoELikelihoodVal(maxMove) > refLikelihoodVal
                   display('Product of experts predictions dominated.');
              end
              
              % This expert has been modified. Let's update the
              % existing data structures to reflect this change.
              display(['PoE likelihood changed from ' num2str(round(refLikelihoodVal)) ' to ' num2str(round(maxNewPoELikelihoodVal(maxMove))) '.']);
              posLikelihoodArr(steps+1) = maxNewPosLikelihoodVal(maxMove) - expertPosLikelihoodVal;
              refLikelihoodVal = maxNewPoELikelihoodVal(maxMove) - gradients(maxMove);
              newMaxExperts = maxExperts{moveItr};
              newMaxNodeCoords = maxNodeCoords{moveItr};
              experts(parseTrees(:,curItr) == expertToMove,:) = newMaxExperts(parseTrees(:,curItr) == expertToMove,:);
              nodeCoords(changedExperts,:) = newMaxNodeCoords(changedExperts,:);
              steps = steps+1;
         end
    end
    
    % Create a gif image out of diagnostic information.
    if steps<maxSteps
         imageArr = imageArr(1:steps+1);
         posLikelihoodArr = posLikelihoodArr(1:steps+1);
         poeLikelihoodArr= poeLikelihoodArr(1:steps+1);
         diffImageArr = diffImageArr(1:steps+1);
    end
    nonemptyIdx = cellfun(@(x) ~isempty(x), diffImageArr);
    diffMax = max(cellfun(@(x) max(max(x)), diffImageArr(nonemptyIdx)));
    if diffMax>0
         diffImageArr(nonemptyIdx) = cellfun(@(x) x/diffMax, diffImageArr(nonemptyIdx), 'UniformOutput', false);
    end

    % Try creating an image and save the image. If image showing
    % fails, we switch back to normal stuff.
    try
         posPadding = (max(posLikelihoodArr) - min(posLikelihoodArr))/4;
         posLimits = [min(posLikelihoodArr) - posPadding, max(posLikelihoodArr) + posPadding];
         if numel(unique(posLimits)) == 1
              posLimits = [posLimits(1)-1, posLimits(1)+1];
         end
         poePadding = (max(poeLikelihoodArr) - min(poeLikelihoodArr))/4;
         poeLimits = [min(poeLikelihoodArr) - poePadding, max(poeLikelihoodArr) + poePadding];
         if numel(unique(poeLimits)) == 1
              poeLimits = [(poeLimits(1)-1), (poeLimits(1) +1)];
         end
         fileName = [folderName '.gif'];
         for stepItr = 1:min(maxSteps, (steps+1))
              figure('Visible', 'off'), hold on;
              axis square
              subplot(2,2,1), imshow(imageArr{stepItr});
              title('Imagined Data')

              subplot(2,2,2)
              title('Likelihood differences')
              if ~isempty(diffImageArr{stepItr})
                   imshow(diffImageArr{stepItr});
              end

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
    catch %#ok<CTCH>
         % Visualization failed possibly due to parallelization. Let's
         % save the data structures for later use.
         save([folderName '.mat'], 'imageArr', 'posLikelihoodArr', 'poeLikelihoodArr', 'diffImageArr');
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
