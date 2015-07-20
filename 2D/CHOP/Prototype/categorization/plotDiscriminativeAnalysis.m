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
function [] = AnalyzePosteriorProbs(datasetName)
    % Create a color array.
    colorArr = {'b', 'g', 'r', 'c', 'm', 'k'}';
    shapeArr = {'.', 'o', 'x', '+', '*', 's', 'd', 'v', '^', '<', '>', 'p', 'h'}';
    plotShapeArr = allcomb(1:numel(colorArr), 1:numel(shapeArr));
    plotShapeArr = plotShapeArr(randperm(size(plotShapeArr,1)), :);
    plotShapeArr = [colorArr(plotShapeArr(:,1)), shapeArr(plotShapeArr(:,2))];
    plotShapeArr = cellfun(@(x,y) [x, y], plotShapeArr(:,1), plotShapeArr(:,2), 'UniformOutput', false);
    
    if ~exist([pwd '/categorization/analysis/' datasetName '/plots'], 'dir')
        mkdir([pwd '/categorization/analysis/' datasetName  '/plots']);
    end

    % Load relevant info.
    load([pwd '/output/' datasetName '/vb.mat']);
    load([pwd '/output/' datasetName '/export.mat']);
    % Extract simple features from images.
    featureDims = cellfun(@(x) numel(x), vocabulary);
    featureDim = sum(featureDims);
    cumSums = int32(cumsum(featureDims));
    cumSums(2:end) = cumSums(1:(end-1));
    cumSums(1) = 0;
    cumSums(end+1) = featureDim;
    allActivations = int32(exportArr(:,[1,4,5]));
    allActivations(:,1) = allActivations(:,1) + cumSums(allActivations(:,2));
    allActivations = unique(allActivations, 'rows');
    allActivations(:,4) = int32(categoryArrIdx(allActivations(:,3)));
    
    numberOfCategories = max(categoryArrIdx);
    categoryArrList = 1:max(categoryArrIdx);
    catDistArr = zeros(featureDim, size(categoryArrList,2));
    % Analyze each feature's discriminative properties.
    for featureItr = 1:featureDim
        activations = allActivations(allActivations(:,1) == featureItr,:);
        catDist = hist(activations(:,4), categoryArrList);
        catDistArr(featureItr,:) = catDist;
    end
    
    %% In the plotting function, we only consider 13 classes.
    if size(catDistArr,2) > size(plotShapeArr,1)
        catDistArr = catDistArr(:,1:size(plotShapeArr,1));
    end
    numberOfCategories = size(catDistArr, 2);
    for nodeItr = 1:size(catDistArr,1)
       catDistArr(nodeItr,:) = catDistArr(nodeItr,:) / sum(catDistArr(nodeItr,:)); 
    end
    
    for levelItr = 1:numel(vocabulary)
        plotArr = cell(numberOfCategories * 3,1);
        levelCatDistArr = catDistArr((cumSums(levelItr) + 1):cumSums(levelItr+1),:);
        plotArr(1:3:size(plotArr,1)) = repmat({(1:size(levelCatDistArr,1))}, numberOfCategories, 1);
        plotArr(2:3:size(plotArr,1)) = mat2cell(levelCatDistArr', ones(numberOfCategories,1), size(levelCatDistArr,1));
        plotArr(3:3:size(plotArr,1)) = plotShapeArr(1:numberOfCategories);
        figure, hold on;
        plot(plotArr{:});
        axis([0, size(levelCatDistArr,1), 0, 1]);
        hold off;
        saveas(gcf, [pwd '/categorization/analysis/' datasetName  '/plots/level' num2str(levelItr) '.png']);
        close all;
    end
    
end