function [] = visualizeJointParts( shownSamples, shownClusters, refinedClusterSamples, refinedClusters, jointPositions, partClusters, numberOfSubParts, partItr, realRFSize, numberOfClusters, uniqueSubParts, newFolder)

  colors = {[1 0 0], [0 1 0], [0 0 1], [0 0 0], [1 0 1], [0 1 1], [0 0.5 1], [1 1 0], [0 1 0.5], [0.5 0 1], [1 0.5 0], [0.5 1 0]};
  pointSymbols = {'o', '*', 'x', '+', 's', 'd', 'v', 'p', '>', '<', '^', '.'};
%% Apply PCA to the vectors (for visualization).
  [~, jointPCA, latent] = princomp(double(shownSamples)); %#ok<PRINCOMP>
  if latent(3) < 0.01
      jointPCA = jointPCA(:,1:2);
  else
      jointPCA = jointPCA(:,1:3);
  end
  figure('Visible', 'off'), hold on;
  maxCount = min(numberOfClusters+1, numel(pointSymbols));
  for clusterItr = 1:min(numberOfClusters+1, numel(pointSymbols))
      relevantSamples = jointPCA(clusterItr == shownClusters, :);
      [uniqueSamples, ~, uniqueSampleCounts] = unique(relevantSamples, 'rows');
      sampleSizes = hist(uniqueSampleCounts, unique(uniqueSampleCounts));
      sampleSizes = round(log(sampleSizes));
      sampleSizes(sampleSizes<1) = 1;

      if clusterItr < maxCount
          color = colors{clusterItr};
      else
          color = [0.5 0.5 0.5];
          continue;
      end

      % Iterate over sizes.
      uniqueSampleSizes = unique(sampleSizes);
      for sizeItr = 1:numel(unique(sampleSizes))
        if size(jointPCA,2) == 3
            plot3(uniqueSamples(sampleSizes == uniqueSampleSizes(sizeItr),1), uniqueSamples(sampleSizes == uniqueSampleSizes(sizeItr),2), uniqueSamples(sampleSizes == uniqueSampleSizes(sizeItr),3), 'LineStyle', 'none', 'Color', color, 'Marker', pointSymbols{clusterItr}, 'MarkerSize', uniqueSampleSizes(sizeItr) * 2);
        else
            plot(uniqueSamples(sampleSizes == uniqueSampleSizes(sizeItr),1), uniqueSamples(sampleSizes == uniqueSampleSizes(sizeItr),2), 'LineStyle', 'none', 'Color', color, 'Marker', pointSymbols{clusterItr}, 'MarkerSize', uniqueSampleSizes(sizeItr) * 2);
        end
      end
  end
  hold off;
  saveas(gcf, [newFolder '/part' num2str(partItr) 'combs.fig']);
  hardcopy(gcf,[newFolder '/part' num2str(partItr) 'combs.eps'],'-depsc2');
  close(gcf);

  %% Print the filters for every part.
  for clusterItr = 1:numberOfClusters
      featureMaps = zeros(numberOfSubParts, realRFSize, realRFSize, 'single');
      discretePos = fix(jointPositions(partClusters == clusterItr,:)) + floor(realRFSize/2) + 1;
      for mapItr = 1:numberOfSubParts
           featureMap = squeeze(featureMaps(mapItr, :, :));
           tmpArr = discretePos(:, ((mapItr-1)*2+1):(mapItr*2));
           tmpIdx = tmpArr(:,1) + (tmpArr(:,2) - 1) * realRFSize;
           [t1,~,t2] = unique(tmpIdx);
           featureMap(t1) = featureMap(t1) + hist(t2, 1:numel(t1))';
           featureMaps(mapItr,:,:) = featureMap;
      end
      featureMaps = double(featureMaps / max(max(max(featureMaps))));
      combinedImg = zeros((realRFSize+2), (1+(numel(uniqueSubParts))*(realRFSize+1)), 3, 'uint8');
      for mapItr = 1:numberOfSubParts
           rgbImg = label2rgb(uint8(255 * squeeze(featureMaps(mapItr,:,:))), 'jet', 'k');
           combinedImg(2:(end-1), ((mapItr-1) * (realRFSize+1) + 2):((mapItr) * (realRFSize+1)), :) = rgbImg;
      end
      combinedImg([1, end], :, :) = 255;
      combinedImg(:, 1:(realRFSize+1):end, :) = 255;
      imwrite(combinedImg, [newFolder '/part' num2str(partItr) '_Variation' num2str(clusterItr) 'All.png']);
  end

  %% Print the filters for every part for refined case.
  for clusterItr = 1:numberOfClusters
      featureMaps = zeros(numberOfSubParts, realRFSize, realRFSize, 'single');
      discretePos = fix(refinedClusterSamples(refinedClusters == clusterItr,:)) + floor(realRFSize/2) + 1;
      for mapItr = 1:numberOfSubParts
           featureMap = squeeze(featureMaps(mapItr, :, :));
           tmpArr = discretePos(:, ((mapItr-1)*2+1):(mapItr*2));
           tmpIdx = tmpArr(:,1) + (tmpArr(:,2) - 1) * realRFSize;
           [t1,~,t2] = unique(tmpIdx);
           featureMap(t1) = featureMap(t1) + hist(t2, 1:numel(t1))';
           featureMaps(mapItr,:,:) = featureMap;
      end
      featureMaps = double(featureMaps / max(max(max(featureMaps))));
      combinedImg = zeros((realRFSize+2), (1+(numel(uniqueSubParts))*(realRFSize+1)), 3, 'uint8');
      for mapItr = 1:numberOfSubParts
           rgbImg = label2rgb(uint8(255 * squeeze(featureMaps(mapItr,:,:))), 'jet', 'k');
           combinedImg(2:(end-1), ((mapItr-1) * (realRFSize+1) + 2):((mapItr) * (realRFSize+1)), :) = rgbImg;
      end
      combinedImg([1, end], :, :) = 255;
      combinedImg(:, 1:(realRFSize+1):end, :) = 255;
      imwrite(combinedImg, [newFolder '/part' num2str(partItr) '_Variation' num2str(clusterItr) 'Refined.png']);
  end

end

