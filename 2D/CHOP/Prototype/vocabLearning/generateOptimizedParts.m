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
function [] = generateOptimizedParts(vocabLevel, vocabulary, levelItr, options)
    vocabulary{levelItr} = vocabLevel;
    filterSize = size(options.filters{1});
    halfSize = ceil(filterSize(1)/2); 
    samplesPerPart = 1;
    batchFlag = true; % Turn on/off batch gradient descent.
    rng('shuffle');
    
    if strcmp(options.filterType, 'gabor')
        inhibitionHalfSize = options.gabor.inhibitionRadius;
        stride = options.gabor.stride;
    else
        inhibitionHalfSize = options.auto.inhibitionRadius;
        stride = options.auto.stride;
    end
    minPixelValue = 1/255;
    filters = options.filters;
    filters = cellfun(@(x) (x - min(min(x))) / (max(max(x)) - min(min(x))), filters, 'UniformOutput', false);
    
    % As a case study, we replace gabors with 1D gaussian filters
    % stretched.
    filterSize = size(filters{1},1);
    vals = normpdf(1:filterSize, (filterSize+1)/2, 5);
    vals = vals/max(vals);
    firstFilter = repmat(vals, filterSize, 1);
    angle = 180/numel(filters);
    for filterItr = 0:(numel(filters)-1)
         curAngle = -angle * filterItr;
         curFilter = imrotate(firstFilter, curAngle, 'bilinear', 'crop');
         curFilter(curFilter<minPixelValue) = minPixelValue;
         filters{filterItr+1} = curFilter;
    end
 %   visFilters = filters;
    
%    H = fspecial('gaussian', 7, 1);
    visFilters = options.filters;
    visFilters = cellfun(@(x) (x - min(min(x))) / (max(max(x)) - min(min(x))), visFilters, 'UniformOutput', false);
    for filterItr = 1:numel(visFilters)
         curFilter = visFilters{filterItr};
 %        curFilter = imfilter(curFilter, H, 'full');
         curFilter(curFilter<minPixelValue) = minPixelValue;
         visFilters{filterItr} = curFilter;
    end
%     
   % Find the correct image size.
   imageSize = options.receptiveFieldSize * stride * (options.poolDim ^ (levelItr-2)) * (inhibitionHalfSize+1) + halfSize * 2;
   prevImageSize = options.receptiveFieldSize * stride * (options.poolDim ^ (levelItr-3)) * (inhibitionHalfSize+1);
   imageSize = [imageSize, imageSize];
   
   % For now, we make the image bigger.
   imageSize = round(imageSize * 1.5);
   
   %% First, for efficiency, we obtain pixel-level predictions for every part
   parfor vocabNodeItr = 1:numel(vocabLevel)
        for sampleItr = 1:samplesPerPart
             % Obtain optimized projections.
             nodes = [vocabNodeItr, 0, 0, levelItr];
             optimizeImagination(nodes, vocabulary, imageSize, prevImageSize, filters, visFilters, sampleItr, batchFlag, options.datasetName);
        end
    end
   
    % Finally, create output images.
    for vocabNodeItr = 1:numel(vocabLevel)
         for sampleItr = 1:samplesPerPart
              folderName = [pwd '/optimization/' options.datasetName '/level' num2str(levelItr) '/' num2str(vocabNodeItr) '_sample' num2str(sampleItr)];
              load([folderName '.mat'], 'imageArr', 'posLikelihoodArr', 'poeLikelihoodArr', 'diffImageArr');
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
              for stepItr = 1:numel(imageArr)
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
                         xlim([1, numel(imageArr)]);
                   end
                   title('Change in position likelihood')

                   subplot(2,2,4), plot(1:stepItr, poeLikelihoodArr(1:stepItr));
                   ylim(poeLimits);
                   if stepItr>1
                         xlim([1, numel(imageArr)]);
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
    end
end