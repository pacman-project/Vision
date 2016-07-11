%> Name: AnalyzePosteriorProbs
%>
%> Description: This function analyzes the posterior probabilities of
%> classes given features. A confusion matrix is created after the process. 
%> For every feature P in the vocabulary, we consider the classes {C_i} where
%> P(C_i|P) > 0. 
%>
%> @param vocabulary Learned vocabulary including category probabilities for 
%> each node.
%> @param exportArr Set of realizations.
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 04.07.2014
function [] = plotDiscriminativeAnalysis(datasetName)
    % Create a color array.
    colorArr = {'b', 'g', 'r', 'c', 'm', 'k'}';
    options = SetParameters(datasetName, 'train');
    shapeArr = {'.', 'o', 'x', '+', '*', 's', 'd', 'v', '^', '<', '>', 'p', 'h'}';
    plotShapeArr = allcomb(1:numel(colorArr), 1:numel(shapeArr));
    [~, plotShapeArrSortIdx] = sort(plotShapeArr(:,2));
    plotShapeArr = plotShapeArr(plotShapeArrSortIdx,:);
    orgPlotShapeArr = plotShapeArr;
    plotShapeArr = [colorArr(plotShapeArr(:,1)), shapeArr(plotShapeArr(:,2))];
    plotShapeArr = cellfun(@(x,y) [x, y], plotShapeArr(:,1), plotShapeArr(:,2), 'UniformOutput', false);
    markerSize = 20;
    
    if ~exist([pwd '/categorization/analysis/' datasetName '/plots'], 'dir')
        mkdir([pwd '/categorization/analysis/' datasetName  '/plots']);
    end
    load([pwd '/output/' datasetName '/vb.mat'], 'categoryNames'); 

    % Load relevant info.
    load([pwd '/categorization/analysis/' datasetName '/discriminativeAnalysis.mat']);
    load([options.currentFolder '/output/' datasetName '/export.mat'], 'categoryArrIdx');
    
    %% PLOT Precision of parts across all layers.
    numberOfCategories = size(precisionArr, 2) - 2;
    numberOfLevels = max(precisionArr(:,1));
    levelPrecisionArr = zeros(1, numberOfLevels);
    for levelItr = 1:numberOfLevels
        plotArr = cell(numberOfCategories * 3,1);
        levelArr = precisionArr(precisionArr(:,1) == levelItr,3:end);
        levelArr(levelArr == 0) = NaN;
        levelPrecisionArr(levelItr) = mean(levelArr(~isnan(levelArr)));
        plotArr(1:3:size(plotArr,1)) = repmat({(1:size(levelArr,1))}, numberOfCategories, 1);
        plotArr(2:3:size(plotArr,1)) = mat2cell(levelArr', ones(numberOfCategories,1), size(levelArr,1));
        plotArr(3:3:size(plotArr,1)) = plotShapeArr(1:numberOfCategories);
        figure, hold on;
  %      plot(plotArr{:}, 'MarkerSize', markerSize);
        
        % Now, we fit 2nd order polynomials to the data for better
        % visualization.
        for categoryItr = 1:numberOfCategories
            x = plotArr{1 + 3 * (categoryItr-1)};
            y = plotArr{2 + 3 * (categoryItr-1)};
            idx = ~isnan(y);
            y = y(idx);
            x = x(idx);
            if numel(x) > 2
                P = polyfit(x, y, 3);
                yfit = polyval(P, x);
                plot(x, yfit, [colorArr{orgPlotShapeArr(categoryItr,1)} '-.']);
            end
        end
        
        axis([0, size(levelArr,1), 0, 100]);
        title('Precision of parts as classifiers.');
        xlabel('Part ids (ordered by MDL scores after inhibition)');
        ylabel('Precision %');
        legend(categoryNames{:});
        hold off;
        saveas(gcf, [pwd '/categorization/analysis/' datasetName  '/plots/level' num2str(levelItr) '_precision.png']);
        close all;
    end
    
    %% PLOT Recall of parts across all layers.
    levelRecallArr = zeros(1, numberOfLevels);
    for levelItr = 1:numberOfLevels
        plotArr = cell(numberOfCategories * 3,1);
        levelArr = recallArr(recallArr(:,1) == levelItr,3:end);
        levelArr(levelArr == 0) = NaN;
        levelRecallArr(levelItr) = mean(levelArr(~isnan(levelArr)));
        plotArr(1:3:size(plotArr,1)) = repmat({(1:size(levelArr,1))}, numberOfCategories, 1);
        plotArr(2:3:size(plotArr,1)) = mat2cell(levelArr', ones(numberOfCategories,1), size(levelArr,1));
        plotArr(3:3:size(plotArr,1)) = plotShapeArr(1:numberOfCategories);
        figure, hold on;
  %      plot(plotArr{:}, 'MarkerSize', markerSize);
        
        % Now, we fit 2nd order polynomials to the data for better
        % visualization.
        for categoryItr = 1:numberOfCategories
            x = plotArr{1 + 3 * (categoryItr-1)};
            y = plotArr{2 + 3 * (categoryItr-1)};
            idx = ~isnan(y);
            y = y(idx);
            x = x(idx);
            if numel(x) > 2
                P = polyfit(x, y, 3);
                yfit = polyval(P, x);
                plot(x, yfit, [colorArr{orgPlotShapeArr(categoryItr,1)} '-.']);
            end
        end
        
        axis([0, size(levelArr,1), 0, 100]);
        title('Recall of parts as classifiers.');
        xlabel('Part ids (ordered by MDL scores after inhibition)');
        ylabel('Recall %');
        legend(categoryNames{:});
        hold off;
        saveas(gcf, [pwd '/categorization/analysis/' datasetName  '/plots/level' num2str(levelItr) '_recall.png']);
        close all;
    end
    
     instanceCounts = hist(categoryArrIdx, unique(categoryArrIdx));
     maxInstanceCountPerClass = max(instanceCounts);
     fscoreBounds = zeros(maxInstanceCountPerClass, 2);
     for itr = 1:maxInstanceCountPerClass
          recCat = itr/maxInstanceCountPerClass;
          preInst = 1/itr;
          fscoreBounds(itr,1) = 2 * recCat / (1 + recCat);
          fscoreBounds(itr,2) = 2 * preInst / (1 + preInst);
     end
     fscoreBounds = 100 * fscoreBounds;
     
    %% PLOT Fscore of parts across all layers.
    levelFscoreArr = zeros(1, numberOfLevels);
    for levelItr = 1:numberOfLevels
        levelArr = fscoreArr(fscoreArr(:,1) == levelItr,3:end);
        levelArr(levelArr == 0) = NaN;
        levelArr(levelArr == 0) = NaN;
        levelFscoreArr(levelItr) = mean(levelArr(~isnan(levelArr)));
        plotArr(1:3:size(plotArr,1)) = repmat({(1:size(levelArr,1))}, numberOfCategories, 1);
        plotArr(2:3:size(plotArr,1)) = mat2cell(levelArr', ones(numberOfCategories,1), size(levelArr,1));
        plotArr(3:3:size(plotArr,1)) = plotShapeArr(1:numberOfCategories);
        figure, hold on;
        
        for categoryItr = 1:numberOfCategories
            x = plotArr{1 + 3 * (categoryItr-1)};
            y = plotArr{2 + 3 * (categoryItr-1)};
            idx = ~isnan(y);
            y = y(idx);
            x = x(idx);
            if numel(x) > 2
                P = polyfit(x, y, 3);
                yfit = polyval(P, x);
                plot(x, yfit, [colorArr{orgPlotShapeArr(categoryItr,1)} '-.']);
            end
        end
 %       partLabels = 1:size(levelArr,1);
        
        % Now, we fit 2nd order polynomials to the data for better
        % visualization.
%         for categoryItr = 1:numberOfCategories
%             validPartIdx = partCategories == categoryItr;
%             if nnz(validPartIdx) == 0
%                  continue;
%             end
%              
%             y = partFscores(validPartIdx);
%             x = partLabels(validPartIdx);
%             z = plotShapeArr{categoryItr};
%              
%             plot(x,y,z, 'MarkerSize', markerSize);
%         end
        
        axis([0, size(levelArr,1), 0, 100]);
        title('Fscore of parts as classifiers.');
        xlabel('Part ids (ordered by MDL scores after inhibition)');
        ylabel('Fscore %');
        legend(categoryNames{:});
        hold off;
        saveas(gcf, [pwd '/categorization/analysis/' datasetName  '/plots/level' num2str(levelItr) '_fscore.png']);
        close all;
    end
    
    %% PLOT Instance based F-scores vs Category F-scores.
   for levelItr = 1:numberOfLevels
        levelFscoreCatlArr = fscoreArr(fscoreArr(:,1) == levelItr,3:end);
        levelFscoreCatlArr = max(levelFscoreCatlArr, [], 2);
        levelFscoreInstArr = fscoreInstanceArr(fscoreInstanceArr(:,1) == levelItr,3) * 100;
        
        figure, hold on;
        plot(levelFscoreInstArr, levelFscoreCatlArr, plotShapeArr{3}, 'MarkerSize', markerSize);
        
        axis([0, 100, 0, 100]);
        plot(fscoreBounds(:,2), fscoreBounds(:,1), 'b-', 'LineWidth', 2);
        title('Instance F-scores vs. Category F-scores');
        ylabel('Category F-scores %');
        xlabel('Instance F-scores %');
        hold off;
        saveas(gcf, [pwd '/categorization/analysis/' datasetName  '/plots/level' num2str(levelItr) '_catInsFscore.png']);
        close all;
   end
   
    %% PLOT Shareability of parts across all layers.
    levelShareabilityArr = zeros(1, numberOfLevels);
    for levelItr = 1:numberOfLevels
        plotArr = cell(3,1);
        levelArr = shareabilityArr(shareabilityArr(:,1) == levelItr,3:end);
        levelArr(levelArr == 0) = NaN;
        levelShareabilityArr(levelItr) = mean(levelArr(~isnan(levelArr)));
        plotArr(1) = {(1:size(levelArr,1))};
        plotArr(2) = {levelArr'};
        plotArr(3) = plotShapeArr(1);
        figure, hold on;
%        plot(plotArr{:}, 'MarkerSize', markerSize);
        
        % Now, we fit 2nd order polynomials to the data for better
        % visualization.
        x = plotArr{1};
        y = plotArr{2};
        idx = ~isnan(y);
        y = y(idx);
        x = x(idx);
        if numel(x) > 2
            P = polyfit(x, y, 3);
            yfit = polyval(P, x);
            plot(x, yfit, [colorArr{orgPlotShapeArr(1)} '-.']);
        end
        
        plot([0, size(levelArr,1)], [100/numberOfCategories, 100/numberOfCategories], 'r-', 'LineWidth', 2);
        axis([0, size(levelArr,1), 0, 100]);
        title('Shareability of parts among categories.');
        xlabel('Part ids (ordered by MDL scores after inhibition)');
        ylabel('Shareability %');
        hold off;
        saveas(gcf, [pwd '/categorization/analysis/' datasetName  '/plots/level' num2str(levelItr) '_shareability.png']);
        close all;
    end
    
    %% Finally, we create a single plot where we show change of average measurements across layers. 
    % Average precision, recall, fscore and shareability of the parts
    % across all layers will be in this plot.
    allDataArr = [levelFscoreArr; levelPrecisionArr; levelRecallArr; levelShareabilityArr];
    plotArr = cell(12,1);
    
    plotArr(1:3:size(plotArr,1)) = repmat({(1:numberOfLevels)}, 4, 1);
    plotArr(2:3:size(plotArr,1)) = mat2cell(allDataArr, ones(4,1), numberOfLevels);
    plotArr(3:3:size(plotArr,1)) = plotShapeArr(1:4);
    figure, hold on;
    plot(plotArr{:}, 'MarkerSize', markerSize);
    
    for itr = 1:4
        x = plotArr{(itr-1)*3 + 1};
        y = plotArr{(itr-1)*3 + 2};
        idx = ~isnan(y);
        y = y(idx);
        x = x(idx);
        if numel(x) > 2
            P = polyfit(x, y, 3);
            yfit = polyval(P, x);
            plot(x, yfit, [colorArr{orgPlotShapeArr(itr)} '-.']);
        end
    end
    axis([1, numberOfLevels, 0, 100]);
    title('Vocabulary properties across layers.');
    xlabel('Layer ids');
    ylabel('Percentage (%)');
    legend({'Fscore', 'Precision', 'Recall', 'Shareability'});

    hold off;
    saveas(gcf, [pwd '/categorization/analysis/' datasetName  '/plots/allMetrics.png']);
    close all;
    
end
