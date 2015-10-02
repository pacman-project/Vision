%> Name: printConvolutionFilters
%>
%> Description: Given the vocabulary nodes in vocabLevel, this function
%> visualizes the convolution filters into the debug folder, as well as
%> saving them into the output folder in .mat format. 
%>
%> @param vocabLevel The vocabulary nodes.
%> @param modes Modes of the previous level.
%> @param levelItr The level id (for visualization).
%> @param debug The debug flag, 1 if debugging.
%> @param Output folder for printing/writing.
%> 
%> @retval none
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 17.09.2015
function [ ] = printConvolutionFilters( vocabLevel, orNodeProbs, modes, modeProbArr, levelItr, previousNodeCount, debug, CNNFolder )
     % We create the folder structure for visualization.
     if debug
          for itr = 1:numel(vocabLevel);
               if ~exist([CNNFolder '/level' num2str(levelItr) '/' num2str(itr)], 'dir')
                    mkdir([CNNFolder '/level' num2str(levelItr) '/' num2str(itr)]);
               end
          end
     end
     
     % Create a distribution image for center nodes.
     convFilterSize = [size(modeProbArr, 2), size(modeProbArr, 3)];
     centerCoords = ceil(convFilterSize/2);
     centerMode = zeros(convFilterSize,  'single');
     centerMode(centerCoords(1), centerCoords(2)) = 1;
     
     % Go through the list of filters and visualize them.
     convFilters = cell(numel(vocabLevel),1);
     parfor vocabLevelItr = 1:numel(vocabLevel)
          convFilter = zeros(previousNodeCount, convFilterSize(1), convFilterSize(2), 'single');
          
          % Add center node's distribution.
          children = vocabLevel(vocabLevelItr).children;
          for childItr = 1:numel(children)
               
               if childItr == 1
                    addedMode = centerMode;
               else
                    relevantModeIdx = find(modes(:,1) == children(1) & modes(:,2) == children(childItr) & modes(:, 3) == vocabLevel(vocabLevelItr).adjInfo(childItr-1,3) , 1, 'first');
                    addedMode =squeeze(modeProbArr(relevantModeIdx,:,:));
               end
               % For every or node, we print the nodes with the given
               % distributions.
               relevantOrNodeProbs = orNodeProbs{children(childItr)};
               for orNodeItr = 1:size(relevantOrNodeProbs,1)
                    tempChild = round(relevantOrNodeProbs(orNodeItr,1));
                    tempMode = addedMode * single(relevantOrNodeProbs(orNodeItr,2));
                    convFilter(tempChild, :, :) = max(squeeze(convFilter(tempChild, :, :)), ...
                         tempMode);
               end
          end
          
          % Finally, we print the distributions.
          printedFolder = [CNNFolder '/level' num2str(levelItr) '/' num2str(vocabLevelItr)];
          convFilters{vocabLevelItr} = convFilter;
          convFilter = round(convFilter * 65535);
          for itr = 1:size(convFilter,1)
               img = squeeze(convFilter(itr,:,:));
               % If it's full of zeros, we move on.
               if nnz(img) == 0
                    continue;
               end
               tempImg = [img, ones(size(img,1),1) * 65535];
               tempImgColored = label2rgb(tempImg, 'jet', 'k');
               img = tempImgColored(:,1:(end-1),:);
               imwrite(img, [printedFolder '/' num2str(itr) '.png']);
          end
     end
     
     % Save convolution filters as well.
     for vocabLevelItr = 1:numel(vocabLevel)
          convFilter = convFilters{vocabLevelItr}; %#ok<NASGU>
          save([CNNFolder '/level' num2str(levelItr) '/' num2str(vocabLevelItr) '/convFilter.mat'], 'convFilter');
     end
end

