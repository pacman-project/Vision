function [ ] = visualizePartStats( partStats, reconstructiveSubs, discriminativeSubs, fscoreSubs, folderName )
     numberOfParts = numel(partStats.fscoreArr);
     numberOfBins = 20;
     names = fieldnames(partStats);
     names = cellfun(@(x) x(1:(end-3)), names, 'UniformOutput', false);
     numberOfFields = numel(names);
     setNames = {'Reconstructive selection', 'MRMR selection', 'F-score selection'};
     
     %% Visualize the data for every field.
     for fieldItr = 1:numberOfFields
          tempArr = getfield(partStats, [names{fieldItr} 'Arr']); %#ok<GFLD>
          
          % Determine histogram min/max values.
          maxVal = max(tempArr);
          halfStep = maxVal/((numberOfBins)*2);
          minVal = halfStep;
          maxVal = maxVal - halfStep;
          histCenters = minVal:(halfStep*2):maxVal;
          limitCount = log2(max(hist(tempArr, histCenters)));
          
          %% Obtain stats for every selection mechanism.
          for subSetItr = 1:3
               % Set correct set.
               if subSetItr == 1
                    selectedSubs = reconstructiveSubs;
               elseif subSetItr == 2
                    selectedSubs = discriminativeSubs;
               else
                    selectedSubs = fscoreSubs;
               end
               eliminatedSubs = setdiff(1:numberOfParts, selectedSubs);
               
               % Collect data for sub set.
               selectedPartStats = tempArr(selectedSubs);
               eliminatedPartStats = tempArr(eliminatedSubs);
               
               % First, we create histogram in log-scale.
               selectedVals = hist(selectedPartStats, histCenters);
               eliminatedVals = hist(eliminatedPartStats,histCenters);
               totalVals = selectedVals + eliminatedVals;
               totalVals = log2(totalVals);
               totalVals(totalVals < 0) = 0;
               selectedVals = log2(selectedVals);
               selectedVals(selectedVals < 0) = 0;
               
               % Next, we visualize the values.
               if subSetItr == 1
                    figure('Visible', 'off'), subplot(1, 3, 1), bar(totalVals, 'r');
               else
                    subplot(1, 3, subSetItr), bar(totalVals, 'r');
               end
               hold on;
               bar(selectedVals, 'g');
               if subSetItr == 2
                    legend('Eliminated Parts', 'Selected Parts');
                    legend boxoff
               end
               ylabel('Number of Parts');
               xlabel(names{fieldItr});
               xlim([0 numberOfBins + 0.5]);
               ylim([0 limitCount]);
               
               % Set axes labels.
               set(gca, 'XTick', 0:floor(numberOfBins/5):numberOfBins);
               labels = num2cell(halfStep * 2 * [0:floor(numberOfBins/5):numberOfBins]);
               if maxVal > 1
                    labels = cellfun(@(x) round(x), labels, 'UniformOutput', false);
                    labels = cellfun(@(x) num2str(x), labels, 'UniformOutput', false);
               else
                    labels = cellfun(@(x) num2str(x, 1), labels, 'UniformOutput', false);
               end
               set(gca, 'XTickLabel', labels');
               tempVals = get(gca, 'YTick');
               tempLabels = cell(numel(tempVals),1);
               for itr = 2:numel(tempVals)
                    tempLabels{itr} = ['2^' num2str(tempVals(itr))];
               end
               set(gca, 'YTickLabel', tempLabels);
               set(gca,'fontsize',6);
               
               % Set title
               title(setNames{subSetItr});
               hold off;
               
              % Make the subplots square.
              axesHandles = findobj(get(gcf,'Children'), 'flat','Type','axes');
              axis(axesHandles,'square');
          end
          suptitle(['Statistics for ' names{fieldItr}]);
          set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
          saveas(gcf, [folderName '/' names{fieldItr} '.png']);
          close gcf;
     end
end

