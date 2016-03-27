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
function [ refModalImg, experts ] = optimizeImagination( nodes, vocabulary, imageSize, prevImageSize, filters, visFilters, sampleItr, ~, datasetName, likelihoodLookupTable, fileName)
     stopVal = 0.01;
%     maxSteps = 10 * size(nodes,1);
     maxSteps = 500;
     minOptimizationLayer = 3;
     minLikelihoodChange = 0.01;
     % If an expert (on average) has better likelihood than this, it means it's more or less agreed.
     likelihoodThr = likelihoodLookupTable(1,1) * 1.005; % Change in likelihood that's enough to reconsider that sub.
     poeCounter = 0;
     
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
     
     % Move parameters.
     movesPerChild = 2;
     maxMoves = 20;
     minMoves = 10;
     
     if isempty(nodes)
          refModalImg = zeros(imageSize, 'uint8');
          experts = [];
          return;
     end
     
     % If poeFlag is true, we are searching for pixel level agreement.
     poeFlag = true;
     % Position flags are the strings that keep things in place.
     positionFlag = false;
     
     % Obtain number of choices per receptive field.
     imageChoices = prevImageSize * prevImageSize;
     
     % Arguments relating to availability of different moves.
     moveFlags = [1; 1; 1] > 0; % 1 for position moves, 2 is for or moves, 3 for rotation moves.
     
     % Shut down warnings.
     warning('off','all');
     
     %% First, imagine layer 1 nodes.
     [ experts, subChildrenExperts, subChildren, orNodeChoices, orNodeChoiceCounts ] = projectNode(nodes, vocabulary, sampleString);
     nodeAngles = zeros(numel(subChildren),1);
     nodeIds = nodes(:,1);
     nodeExpertCounts = cellfun(@(x) sum(cellfun(@(y) size(y,1), x)), subChildrenExperts);
     nodeExpertIds = zeros(size(experts,1),1);
     startOffset = 1;
     for tempItr = 1:size(nodes,1)
          nodeExpertIds(startOffset:(startOffset + nodeExpertCounts(tempItr) - 1)) = tempItr;
          startOffset = startOffset + nodeExpertCounts(tempItr);
     end
     experts = double(experts);
     
    % If rotation flag is on, we'll consider a wider range of filters.
    numberOfFilters = size(filters,3);
    rotationFlag = numberOfFilters > numel(vocabulary{1});
    numberOfRealFilters = numel(vocabulary{1});
    if numberOfFilters > numel(vocabulary{1})
          filterIds = round(((180/numberOfRealFilters) * (0:(numberOfRealFilters-1))) / (180/numberOfFilters))' + 1;
          experts(:,1) = filterIds(experts(:,1));
    end
     
    %% Second, we start a gradient-descent procedure to move things around so they match.
    modalImg = [];
    likelihoodMat = [];
     
    % Only move nodes from level below.
    steps = 0;
    curLevelItr = nodes(1, 4);
    
    % TODO REMOVE.
    minOptimizationLayer = curLevelItr;
    imageArr = cell(maxSteps,1);
    posLikelihoodArr = zeros(maxSteps,1);
    poeLikelihoodArr = zeros(maxSteps,1);
    
    %% Continue with gradient descent until optimized.
    while steps < maxSteps && curLevelItr >= minOptimizationLayer
         % Get appropriate vocabulary level.
         currentLevel = vocabulary{curLevelItr};
         numberOfHighExperts = numel(subChildren);
         
         % Select experts to move.
         modifiedExperts = 1:numberOfHighExperts;
         
         %% We obtain the reference modal image and likelihood mat, in order to be able to select parts based on their local likelihood values.
         [refModalImg, ~, ~] = obtainPoE(experts, modalImg, likelihoodMat, imageSize, visFilters, []);
         imageArr(steps+1) = {refModalImg};
         [refTrueModalImg, refLikelihoodMat, refLikelihoodVal] = obtainPoE(experts, modalImg, likelihoodMat, imageSize, filters, likelihoodLookupTable);
         poeCounter = poeCounter  + 1;
         
         % Check which experts are to be moved. We'll move the expert which
         % has lowest average likelihood.
         avgLikelihoods = zeros(numel(modifiedExperts),1);
         expertIdx = cell(numel(modifiedExperts),1);
         zeroMatrix = zeros(imageSize) > 0;
         for expertItr = 1:numel(modifiedExperts)
              subExperts = subChildrenExperts{expertItr};
              subExperts = cat(1, subExperts{:});
              lowLevelExperts = sub2ind(imageSize, subExperts(:,2), subExperts(:,3));
              tempMatrix = zeroMatrix;
              tempMatrix(lowLevelExperts) = 1;
              tempOverlapMatrix = imdilate(tempMatrix, c_mask);
              expertIdx(expertItr) = {find(tempOverlapMatrix)};
              avgLikelihoods(expertItr) = mean(refLikelihoodMat(tempOverlapMatrix));
         end
         moveFlagArr = avgLikelihoods < likelihoodThr;
         
         % Find which expert to move.
         [~, sortedExpertIdx] = sort(avgLikelihoods);
         
         %% Iterate over list of experts and move them around.
         while steps<maxSteps
              % If no experts could be moved, move on.
              if nnz(moveFlagArr) == 0
                   break;
              end
              
              % Select the expert with minimum average likelihood.
              curExpertItr = find(moveFlagArr(sortedExpertIdx), 1, 'first');
              tempIdx = sortedExpertIdx(curExpertItr);
              
              % Obtain expert information.
              expertToMove = modifiedExperts(tempIdx);
              expertLabelId = nodeIds(expertToMove);
              expertNode = currentLevel(expertLabelId);
              
              % Get OR node info.
              expertOrNodeChoice = orNodeChoices(expertToMove);
              expertOrNodeChoiceCount = orNodeChoiceCounts(expertToMove);
              
               % First of all, we select the relevant combination of nodes.
               expertChildren = subChildren{expertToMove};
               
               % Calculate pos probability denominator.
               posProbDenom  = imageChoices ^ (size(expertChildren,1) - 1);

               % Then, we obtain the position pdfs.
               childrenPosDistributions = expertNode.childrenPosDistributions{1};
               if numel(childrenPosDistributions) > 1
                    childrenPosDistributions = childrenPosDistributions{expertOrNodeChoice}; %#ok<FNDSB>
               else
                    childrenPosDistributions = childrenPosDistributions{1}; %#ok<FNDSB>
               end
               
              %% We sample moves using available move types.
              stochasticMoves = round(max(minMoves, min(movesPerChild^size(expertChildren,1), maxMoves)));
              moves = generateMoves(stochasticMoves, numel(expertNode.children), moveFlags, expertOrNodeChoice, expertOrNodeChoiceCount);
              numberOfMoves = size(moves,1);
               
%               %% Calculate existing prediction's position likelihood.
%               childExpertCenter = expertChildren(1);
%               if numel(expertChildren) > 1
%                    childExpertPers = expertChildren(2:end);
%               else
%                    childExpertPers = [];
%               end
              % If positions are involved, 
              if positionFlag
                   expertPosLikelihoodVal = GetTreeLikelihood(childExpertCenter, childExpertPers, nodeCoords, childrenPosDistributions, posProbDenom);
              else
                   expertPosLikelihoodVal = 0;
              end
              
              %% Update expert positions based on the perturbations.
              gradients = zeros(numberOfMoves,1);
%              parentExpertCenter = nodeCoords(expertToMove, :);
              
              % Before we start moving, we calculate an initial product of
              % experts by not taking the expert in question into account.
              if poeFlag              
                   removedExperts = -1 * (nodeExpertIds == expertToMove);
                   if nnz(removedExperts == 0) > 0
                        [preModalImg, preLikelihoodMat, ~] = obtainPoE(experts, refTrueModalImg, refLikelihoodMat, imageSize, filters, likelihoodLookupTable, removedExperts);
                        poeCounter = poeCounter + 1;
                   else
                        preModalImg = [];
                        preLikelihoodMat = [];
                   end
              end
              
              % Allocate space for max experts.
              maxExperts = cell(numberOfMoves,1);
              maxSubChildrenExperts = cell(numberOfMoves,1);
              maxSubChildren = cell(numberOfMoves,1);
              maxOrNodeChoices = zeros(numberOfMoves,1);
              maxAngles = zeros(numberOfMoves,1);
              maxNodeExpertIds = cell(numberOfMoves,1);
              maxLikelihoodMat = cell(numberOfMoves,1);
              maxNewPoELikelihoodVal = zeros(numberOfMoves,1);
              maxNewPosLikelihoodVal = zeros(numberOfMoves,1);
              
              % Calculate moves, and their gradients.              
              for moveItr = 1:numberOfMoves
                   if moves(moveItr,1) == 2
                         newNodes = nodes(expertToMove,:);
                         [ newExperts, newSubChildrenExperts, newSubChildren, ~, ~ ] = projectNode(newNodes, vocabulary, sampleString, moves(moveItr,2));
                         newSubChildrenExperts = newSubChildrenExperts{1};
                         newSubChildren = newSubChildren{1};
                   else
                         newExperts = subChildrenExperts{expertToMove};
                         newExperts = cat(1, newExperts{:});
                         newSubChildrenExperts = subChildrenExperts{expertToMove};
                         newSubChildren = subChildren{expertToMove};
                   end
                   if rotationFlag 
                         filterIds = round(((180/numberOfRealFilters) * (0:(numberOfRealFilters-1))) / (180/numberOfFilters))' + 1;
                         newExperts(:,1) = filterIds(newExperts(:,1));
                   end
                   
                   
                   %% Apply the moves, and get a new set of experts.
                   % Or node choice don't forget!
                   [ newExperts, newSubChildrenExperts, newSubChildren, newAngle ] = ...
                        applyMove(nodes(expertToMove, 2:3), newExperts, newSubChildren, newSubChildrenExperts,...
                        moves(moveItr,:), nodeAngles(expertToMove), numberOfRealFilters, numberOfFilters);
                   
                   % Insert new experts to their appropriate places.
                   nodeExpertCounts = cellfun(@(x) sum(cellfun(@(y) size(y,1), x)), subChildrenExperts);
                   nodeExpertIds = zeros(size(experts,1),1);
                   startOffset = 1;
                   for tempItr = 1:size(nodes,1)
                        nodeExpertIds(startOffset:(startOffset + nodeExpertCounts(tempItr) - 1)) = tempItr;
                        startOffset = startOffset + nodeExpertCounts(tempItr);
                   end
                   beforeIdx = nodeExpertIds < expertToMove;
                   beforeExperts = experts(beforeIdx, :);
                   afterIdx = nodeExpertIds > expertToMove;
                   afterExperts = experts(afterIdx, :);
                   numberOfNewExperts = size(newExperts,1);
                   newExperts = cat(1, beforeExperts, newExperts, afterExperts);
                   
                   % Update expert-node associations.
                   addedExperts = zeros(size(newExperts,1),1);
                   addedExperts((size(beforeExperts,1) + 1):(size(beforeExperts,1) + numberOfNewExperts)) = 1;
                   newNodeExpertIds = cat(1, nodeExpertIds(beforeIdx), ones(numberOfNewExperts,1) * expertToMove, nodeExpertIds(afterIdx));
                   
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
                             [~, newLikelihoodMat, ~] = obtainPoE(newExperts, preModalImg, preLikelihoodMat, imageSize, filters, likelihoodLookupTable, addedExperts);
                             poeCounter = poeCounter + 1;
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
                        maxExperts(moveItr) = {newExperts};
                        maxSubChildrenExperts(moveItr) = {newSubChildrenExperts};
                        maxSubChildren(moveItr) = {newSubChildren};
                        maxOrNodeChoices(moveItr) = moves(moveItr,2);
                        maxNodeExpertIds(moveItr) = {newNodeExpertIds};
                        maxAngles(moveItr) = newAngle;
                        maxLikelihoodMat(moveItr) = {newLikelihoodMat};
                        maxNewPoELikelihoodVal(moveItr) = newPoELikelihoodVal;
                        maxNewPosLikelihoodVal(moveItr) = newPosLikelihoodVal;
                        gradients(moveItr) = newDiffVal;
                        break;
                   end
              end
              
              if nnz(gradients) == 0
                   % This expert didn't move, let's move on to the next
                   % one.
                   moveFlagArr(expertToMove) = 0;
                   continue;
              end
              
              % We have moved the expert, now let's update the data
              % structures.
              [~, maxMove] = max(gradients);
              
              % If we're not performing batch gradient descent, let's
              % break and update our reference likelihood values. 
              display(['Expert ' num2str(expertToMove) ' has been modified.']);
              display(['Position likelihood changed from ' num2str(round(expertPosLikelihoodVal)) ' to ' num2str(round(maxNewPosLikelihoodVal(maxMove))) '.']);
              if maxNewPosLikelihoodVal(maxMove) > expertPosLikelihoodVal && maxNewPoELikelihoodVal(maxMove) > refLikelihoodVal
                   display('Both likelihoods improved!');
              elseif maxNewPosLikelihoodVal(maxMove) > expertPosLikelihoodVal
                   display('Position distributions pulled the strings.');
              elseif maxNewPoELikelihoodVal(maxMove) > refLikelihoodVal
                   display('Product of experts predictions dominated.');
              end
              display(['PoE likelihood changed from ' num2str(round(refLikelihoodVal)) ' to ' num2str(round(maxNewPoELikelihoodVal(maxMove))) '.']);
              
              % This expert has been modified. Let's update the
              % existing data structures to reflect this change.
              posLikelihoodArr(steps+1) = maxNewPosLikelihoodVal(maxMove) - expertPosLikelihoodVal;
              poeLikelihoodArr(steps+1) = maxNewPoELikelihoodVal(maxMove);
              refLikelihoodVal = maxNewPoELikelihoodVal(maxMove);
              
              %% Assign max data structures.
              prevExperts = experts;
              experts = maxExperts{maxMove};
              subChildrenExperts{expertToMove} = maxSubChildrenExperts{maxMove};
              subChildren{expertToMove} = maxSubChildren{maxMove};
              orNodeChoices(expertToMove) = maxOrNodeChoices(maxMove);
              nodeExpertIds = maxNodeExpertIds{maxMove};
              nodeAngles(expertToMove) = maxAngles(maxMove);
              steps = steps+1;
              
              %% Calculate new reference image.
              expertTypes = [removedExperts; ones(nnz(addedExperts),1)];
              updatedExperts = [prevExperts; experts(addedExperts>0,:)];
              imageArr(steps+1) = {refModalImg};         
              [refModalImg, ~, ~] = obtainPoE(updatedExperts, refModalImg, [], imageSize, visFilters, [], expertTypes);
              poeCounter = poeCounter + 1;
              refLikelihoodMat = maxLikelihoodMat{maxMove};
              
              %% Update the likelihood ordering array and linear indices.
              tempMatrix = zeroMatrix;
              lowLevelExperts = sub2ind(imageSize, experts(nodeExpertIds == maxMove, 1), experts(nodeExpertIds == maxMove, 2));
              tempMatrix(lowLevelExperts) = 1;
              tempOverlapMatrix = imdilate(tempMatrix, c_mask);
              expertIdx(expertToMove) = {find(tempOverlapMatrix)};
              
              % Update data structures
              oldAvgLikelihoods = avgLikelihoods;
              for expertItr = 1:numel(modifiedExperts)
                   avgLikelihoods(expertItr) = mean(refLikelihoodMat(expertIdx{expertItr}));
              end
              
              % Check which experts' values have been changed.
              moveFlagArr = avgLikelihoods < likelihoodThr | abs(oldAvgLikelihoods-avgLikelihoods) >= minLikelihoodChange;
              
              % Sort likelihoods again.
              [~, sortedExpertIdx] = sort(avgLikelihoods);
         end
         
         % Finished with this level, let's move on to the next one.
         curLevelItr = curLevelItr - 1;
    end
    
    % Create a gif image out of diagnostic information.
    if steps<maxSteps
         imageArr = imageArr(1:steps);
         posLikelihoodArr = posLikelihoodArr(1:steps);
         poeLikelihoodArr= poeLikelihoodArr(1:steps);
    end

    % Try creating an image and save the image. If image showing
    % fails, we switch back to normal stuff.
    display('Writing output!');
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
              subplot(1,2,1), imshow(imageArr{stepItr});
              set(gca,'Position',[0.05 0.05 0.5 0.9])
               title('Imagined Data')
% 
%               subplot(1,2,2), plot(1:stepItr, posLikelihoodArr(1:stepItr));
%               ylim(posLimits);
%               if stepItr>1
%                     xlim([1, min(maxSteps, (steps+1))]);
%               end
%               title('Change in position likelihood')
     %         set(gca,'Position',[0.1 .1 0.75 0.85])
              subplot(1,2,2), plot(1:stepItr, poeLikelihoodArr(1:stepItr));
              set(gca,'Position',[0.62 0.05 0.34 0.8])
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
         save([folderName '.mat'], 'imageArr', 'posLikelihoodArr', 'poeLikelihoodArr');
    end
    
    % Turn warnings on.
     warning('on','all');
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