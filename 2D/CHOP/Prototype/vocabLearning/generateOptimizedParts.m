%> Name: postProcessParts
%>
%> Description: Finds the distance between the parts in vocabLevel.
%> Additionally, it assumes that we have already removed/inhibited some nodes
%> from graphLevel (object graphs), so it removes parts which do not have any
%> instances from vocabLevel. Once they are removed, the parts are ordered by
%> the number of occurences, while the original ordering is preserved in
%> graphLabelAssgnArr to be used in inference. 
%>
%> @param vocabLevel Vocabulary level to be processed.
%> @param graphLevel Object graphs, encoded in a node + adjacency list
%> fashion.
%> @param nodeDistanceMatrix The distance matrix of previous layer's nodes.
%> @param options Program options.
%>
%> @retval vocabLevel Processed vocabulary level.
%> @retval graphLevel Remaining nodes for object graphs, encoded in a 
%> node + adjacency list fashion.
%> @retval newDistanceMatrix Generated distance matrix of this layer's nodes.
%> @retval graphLabelAssgnArr Original ordering of parts in vocabLevel 
%> before occurence-based reordering.
%>
%> @retval lowestCost Minimum matching score.
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 06.05.2014
%> Update on 23.02.2015 Added comments, performance boost.
%> Update on 25.02.2015 Added support for single node subs.
function [] = generateOptimizedParts(vocabulary, levelItr, options)
    vocabLevel = vocabulary{levelItr};
    samplesPerPart = 5;
    rng('shuffle');
    
    % Create data structures required for optimization.
    if ~exist([pwd '/filters/optimizationFilters.mat'], 'file')
         [rfSizes, visFilters, optimizedFilters, likelihoodLookupTable] = createOptimizationStructures(options, levelItr, true);
         save([pwd '/filters/optimizationFilters.mat'], 'rfSizes', 'visFilters', 'optimizedFilters', 'likelihoodLookupTable');
    else
         load([pwd '/filters/optimizationFilters.mat'], 'rfSizes', 'visFilters', 'optimizedFilters', 'likelihoodLookupTable');
    end
    
   % Find the correct image size.
   imageSize = getRFSize(options, levelItr);
   
   % For now, we make the image bigger.
   imageSize = round(imageSize * 1.2);
   vocabulary = vocabulary(1:levelItr);
   datasetName = options.datasetName;
   
    %% Create optimization options.
     optimizationOptions.stopVal = 0.1;
     optimizationOptions.maxSteps = 20;
     optimizationOptions.minOptimizationLayer = 3;
     optimizationOptions.minLikelihoodChange = 0.000001;
     optimizationOptions.movesPerChild = 3;
     optimizationOptions.maxMoves = 10;
     optimizationOptions.minMoves = 5;
     optimizationOptions.likelihoodChangeThr =  1.000001;
     optimizationOptions.moveFlags = [1,1,0]; % 1 for position moves, 2 is for or moves, 3 for rotation moves.
     % If poeFlag is true, we are searching for pixel level agreement.
     optimizationOptions.poeFlag = true;
     % Position flags are the strings that keep things in place.
     optimizationOptions.positionFlag = true;
   
   %% First, for efficiency, we obtain pixel-level predictions for every part
 %  for vocabNodeItr = 61:numel(vocabLevel)
  % for vocabNodeItr = [6, 8, 15, 36, 59, 81, 87, 90, 98, 121]
   printedNodes = [78];
   for vocabNodeItr = 1:numel(printedNodes)
%   for vocabNodeItr = 90
        printedVocabNode = printedNodes(vocabNodeItr);
        for sampleItr = 1:samplesPerPart
             % Obtain optimized projections.
             nodes = [printedVocabNode, round(imageSize(1)/2), round(imageSize(2)/2), levelItr];
             optimizeImagination(nodes, vocabulary, imageSize, rfSizes, optimizedFilters, visFilters, sampleItr, datasetName, likelihoodLookupTable, options, optimizationOptions);
        end
    end
   
%     % Finally, create output images.
%     for vocabNodeItr = 1:numel(vocabLevel)
%          for sampleItr = 1:samplesPerPart
%               folderName = [pwd '/optimization/' options.datasetName '/level' num2str(levelItr) '/' num2str(vocabNodeItr) '_sample' num2str(sampleItr)];
%               load([folderName '.mat'], 'imageArr', 'posLikelihoodArr', 'poeLikelihoodArr', 'diffImageArr');
%               posPadding = (max(posLikelihoodArr) - min(posLikelihoodArr))/4;
%               posLimits = [min(posLikelihoodArr) - posPadding, max(posLikelihoodArr) + posPadding];
%               if numel(unique(posLimits)) == 1
%                    posLimits = [posLimits(1)-1, posLimits(1)+1];
%               end
%               poePadding = (max(poeLikelihoodArr) - min(poeLikelihoodArr))/4;
%               poeLimits = [min(poeLikelihoodArr) - poePadding, max(poeLikelihoodArr) + poePadding];
%               if numel(unique(poeLimits)) == 1
%                    poeLimits = [(poeLimits(1)-1), (poeLimits(1) +1)];
%               end
%               fileName = [folderName '.gif'];
%               for stepItr = 1:numel(imageArr)
%                    figure('Visible', 'off'), hold on;
%                    axis square
%                    subplot(2,2,1), imshow(imageArr{stepItr});
%                    title('Imagined Data')
% 
%                    subplot(2,2,2)
%                    title('Likelihood differences')
%                    if ~isempty(diffImageArr{stepItr})
%                         imshow(diffImageArr{stepItr});
%                    end
% 
%                    subplot(2,2,3), plot(1:stepItr, posLikelihoodArr(1:stepItr));
%                    ylim(posLimits);
%                    if stepItr>1
%                          xlim([1, numel(imageArr)]);
%                    end
%                    title('Change in position likelihood')
% 
%                    subplot(2,2,4), plot(1:stepItr, poeLikelihoodArr(1:stepItr));
%                    ylim(poeLimits);
%                    if stepItr>1
%                          xlim([1, numel(imageArr)]);
%                    end
%                    title('Product of experts likelihood')
% 
%                    hold off;
%                    saveas(gcf,  [folderName '_temp.png']);
%                    im=imread([folderName '_temp.png']);
%                    [imind,cm] = rgb2ind(im, 256);
% 
%                    if stepItr == 1
%                         imwrite(imind, cm, fileName, 'LoopCount', inf, 'DelayTime',2);
%                    else
%                         imwrite(imind, cm, fileName, 'WriteMode', 'append');
%                    end
%                    close(gcf);
%               end
%          end
%     end
end