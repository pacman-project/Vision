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
     stopVal = 1;
     maxSteps = 5 * size(nodes,1);
     minLikelihoodChange = 0.002;
     likelihoodThr = likelihoodLookupTable(1,1) * 1.005; 
     poeCounter = 0;
     % If an expert (on average) has better likelihood than this, it means it's more or less agreed.
     
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
     [experts, parseTrees, nodeIds, nodeCoords, orNodeChoices, orNodeChoiceCounts] = projectNode(nodes, vocabulary, sampleString);
     nodeAngles = zeros(numel(nodeIds),1);
     experts = double(experts);
     
     % To work with more precise coordinates, we use imprecise coordinates
     % here.
     experts(:,2:3) = nodeCoords(unique(parseTrees(:,end)),:);
    
    % If rotation flag is on, we'll consider a wider range of filters.
    numberOfFilters = size(filters,3);
    numberOfRealFilters = numel(vocabulary{1});
    if numberOfFilters > numel(vocabulary{1})
          filterIds = round(((180/numberOfRealFilters) * (0:(numberOfRealFilters-1))) / (180/numberOfFilters))' + 1;
          experts(:,1) = filterIds(experts(:,1));
    end
     
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
    likelihoodMat = [];
     
    % Only move nodes from level below.
    steps = 0;
    curItr = 1;
    imageArr = cell(maxSteps,1);
    posLikelihoodArr = zeros(maxSteps,1);
    poeLikelihoodArr = zeros(maxSteps,1);
    
    %% Continue with gradient descent until optimized.
    while steps < maxSteps && curItr <= size(parseTrees,2)-1
         % Get appropriate vocabulary level.
         curLevelItr = (size(parseTrees,2) - curItr) + 1;
         currentLevel = vocabulary{curLevelItr};
         
         % Select experts to move.
         [modifiedExperts, ~] = unique(parseTrees(:,curItr), 'first');
         modifiedExperts = modifiedExperts';
         
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
         ind = sub2ind(imageSize, experts(:,2), experts(:,3));
         for expertItr = 1:numel(modifiedExperts)
              tempMatrix = zeroMatrix;
              lowLevelExperts = ind(parseTrees(:,curItr) == modifiedExperts(expertItr));
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
              
              % Move the expert!
              expertToMove = modifiedExperts(tempIdx);
              expertLabelId = nodeIds(expertToMove);
              expertNode = currentLevel(expertLabelId);
              
              % Get OR node info.
              expertOrNodeChoice = orNodeChoices(expertToMove);
              expertOrNodeChoiceCount = orNodeChoiceCounts(expertToMove);
              
               % First of all, we select the relevant combination of nodes.
               expertChildren = unique(parseTrees(parseTrees(:, curItr) == expertToMove, curItr+1), 'first');
               childrenLabelCombination = nodeIds(expertChildren);
               relevantCombination = find(ismember(childrenLabelCombination, expertNode.childrenLabelDistributions(:,1:(end-1)), 'rows'));
               
               % Calculate pos probability denominator.
               posProbDenom  = imageChoices ^ (numel(expertChildren));

               % Then, we obtain the position pdfs.
               childrenPosDistributions = expertNode.childrenPosDistributions{1}{relevantCombination}; %#ok<FNDSB>
               
              %% We sample moves using available move types.
              stochasticMoves = round(max(minMoves, min(movesPerChild^numel(expertChildren), maxMoves)));
              moves = generateMoves(stochasticMoves, numel(expertNode.children), moveFlags, expertOrNodeChoice, expertOrNodeChoiceCount);
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
                   removedExperts = -double((parseTrees(:, curItr) == expertToMove));
                   if nnz(removedExperts == 0) > 0
                        [preModalImg, preLikelihoodMat, ~] = obtainPoE(experts, refTrueModalImg, refLikelihoodMat, imageSize, filters, likelihoodLookupTable, removedExperts);
                        poeCounter = poeCounter + 1;
                   else
                        preModalImg = [];
                        preLikelihoodMat = [];
                   end
              end
              addedExperts = double(parseTrees(:, curItr) == expertToMove);
              
              % Allocate space for max experts.
              maxExperts = cell(numberOfMoves,1);
              maxLikelihoodMat = cell(numberOfMoves,1);
              maxNewPoELikelihoodVal = zeros(numberOfMoves,1);
              maxNewPosLikelihoodVal = zeros(numberOfMoves,1);
              maxNodeCoords = cell(numberOfMoves,1);
              
              % Calculate moves, and their gradients.              
              for moveItr = 1:numberOfMoves
                   %% Apply the moves, and get a new set of experts.
                   [tempExperts, tempParseTrees, tempNodeIds, tempNodeCoords, tempOrNodeChoices, tempOrNodeChoiceCounts] = ...
                        applyMove(expertToMove, parentExpertCenter, expertChildren, experts, parseTrees, nodeIds, ...
                        nodeCoords, moves(moveItr,:), nodeAngles, orNodeChoices, orNodeChoiceCounts, nodes, vocabulary, sampleString, numberOfFilters, curItr);
                  
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
                        maxExperts(moveItr) = {tempExperts};
                        maxLikelihoodMat(moveItr) = {newLikelihoodMat};
                        maxNewPoELikelihoodVal(moveItr) = newPoELikelihoodVal;
                        maxNewPosLikelihoodVal(moveItr) = newPosLikelihoodVal;
                        maxNodeCoords(moveItr) = {tempNodeCoords};
                        gradients(moveItr) = newDiffVal;
                        break;
                   end
              end
              
              if nnz(gradients) == 0
                   % This expert didn't move, let's move on to the next
                   % one.
                   moveFlagArr(tempIdx) = 0;
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
              poeLikelihoodArr(steps+1) = maxNewPoELikelihoodVal(maxMove);
              refLikelihoodVal = maxNewPoELikelihoodVal(maxMove);
              newMaxExperts = maxExperts{maxMove};
              newMaxNodeCoords = maxNodeCoords{maxMove};
              prevExperts = experts;
              experts(parseTrees(:,curItr) == expertToMove,:) = newMaxExperts(parseTrees(:,curItr) == expertToMove,:);
              nodeCoords(changedExperts,:) = newMaxNodeCoords(changedExperts,:);
              steps = steps+1;
              
              %% Calculate new reference image.
              % TODO: Change
              expertTypes = [removedExperts; ones(nnz(addedExperts),1)];
              updatedExperts = [prevExperts; experts(addedExperts>0,:)];
              imageArr(steps+1) = {refModalImg};         
              [refModalImg, ~, ~] = obtainPoE(updatedExperts, refModalImg, [], imageSize, visFilters, [], expertTypes);
              poeCounter = poeCounter + 1;
              refLikelihoodMat = maxLikelihoodMat{maxMove};
              
              %% Update the likelihood ordering array and linear indices.
              tempMatrix = zeroMatrix;
              lowLevelExperts = ind(parseTrees(:,curItr) == expertToMove);
              tempMatrix(lowLevelExperts) = 1;
              tempOverlapMatrix = imdilate(tempMatrix, c_mask);
              expertIdx(tempIdx) = {find(tempOverlapMatrix)};
              
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
         curItr = curItr + 1;
    end
    
    % Create a gif image out of diagnostic information.
    if steps<maxSteps
         imageArr = imageArr(1:steps+1);
         posLikelihoodArr = posLikelihoodArr(1:steps+1);
         poeLikelihoodArr= poeLikelihoodArr(1:steps+1);
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
