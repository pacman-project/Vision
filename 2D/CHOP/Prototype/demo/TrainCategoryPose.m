%> Name: TrainCategoryPose
%>
%> Description: Trains the category and pose models for a given dataset.
%>
%> @param datasetName Name of the dataset.
%> 
%> @param exportArr The array consisting of all realizations in mainGraph, of
%> the form: detections Detection array of the form:
%> [labelId1, x1, y1, levelId1, imageId1;
%> [labelId2, x2, y2, levelId2, imageId2;
%> ...]
%> 
%> @param categoryArr nx1 cell array (n = number of Images). Each row
%> corresponds to the category label of the corresponding image, indexed with
%> imageId from exportArr.
%> 
%> @param poseArr nx1 cell array (n = number of Images). Each row
%> corresponds to the pose label of the corresponding image, indexed with
%> imageId from exportArr.
%>
%> @retval models Learned category/pose models to be used in predictions.
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 14.05.2014
function [ models ] = TrainCategoryPose( datasetName, exportArr, categoryArr, poseArr )
    %% Set learning parameters.
    feature_params.binSize=32;
    feature_params.cellSize=32;
	feature_params.numOrient=32;
    feature_params.histSize=64;
    feature_params.feature_type='all';
    feature_params.integration_levels=[1,2,3,4,5];
    feature_params.learnedLevel = max(feature_params.integration_levels);
    feature_params.isTesting = 0;
    feature_params.testImageSize = [];
    feature_params.categoryLevel = 4;
    
    % Create an index-based category array.
    imageIds = exportArr(:,5);
    categoryArrIdx = zeros(numel(categoryArr), 1);
    categoryStrs = {'Mug', 'Bowl', 'Tool', 'Pan'};
    categorySamples = zeros(4,1);
    for categoryItr = 1:4
       samples = cellfun(@(x) strcmp(x, categoryStrs{categoryItr}), categoryArr);
       categoryArrIdx(samples) = categoryItr;
       categorySamples(categoryItr) = nnz(samples);
    end
    feature_params.imageSizes = {[180 220], [180 220], [285 360], [285 360]};
    feature_params.categorySamples = categorySamples;
    
    %% Learn size distributions for every class in order to eliminate irrelevant windows in testing.
    minCounts = zeros(4,1);
    meanCounts = zeros(4,1);
    for categoryItr = 1:4
        samples = find(categoryArrIdx==categoryItr);
        numberOfSamples = numel(samples);
        realizationCounts = zeros(numberOfSamples,1);
        for sampleItr = 1:numberOfSamples
            realizationCounts(sampleItr) = nnz(imageIds == samples(sampleItr));
        end
        minCounts(categoryItr) = round(min(realizationCounts));
        meanCounts(categoryItr) = round(mean(realizationCounts));
    end
    feature_params.minCounts = minCounts;
    feature_params.meanCounts = meanCounts;
    
    %% Learn contributions of high-level compositions to categories.
    avgRelativePoints = zeros(4,2);
    for categoryItr = 1:4
        samples = find(categoryArrIdx==categoryItr);
        avgPoints = 0;
        avgPointCtr = 0;
        for sampleItr = 1:numel(samples)
            reals = exportArr(exportArr(:,5) == samples(sampleItr) & exportArr(:,4) == 1, 2:3);
            if ~isempty(reals)
                avgPoint = mean(reals, 1);
                avgPoints = avgPoints + avgPoint;
                avgPointCtr = avgPointCtr + 1;
            end
        end
        avgRelativePoints(categoryItr, :) = (avgPoints / avgPointCtr) - (feature_params.imageSizes{categoryItr}/2);
    end
    feature_params.avgRelativePoints = avgRelativePoints;
    
    %% Find contribution of each part to the categories.
    % Learning category level nodes.
%     categoryLevels = unique(exportArr(:,4));
%     categoryLevels = categoryLevels(categoryLevels>=feature_params.categoryLevel)';
%     contributions = cell(numel(categoryLevels),1);
%     for levelItr = 1:numel(categoryLevels)
%         realizations = exportArr(exportArr(:,4) == categoryLevels(levelItr),:);
%         partIds = unique(realizations(:,1));
%         
%         partContributions = zeros(numel(partIds), 6);
%         for partItr = 1:numel(partIds)
%             imageIds = realizations(realizations(:,1) == partIds(partItr), 5);
%             imageCategories = categoryArrIdx(imageIds);
%             probs = hist(imageCategories, 1:4);
%             probs = probs / sum(probs);
%             partContributions(partItr,:) = [partIds(partItr), categoryLevels(levelItr), probs/numel(partIds)];
%         end
%         contributions(levelItr) = {partContributions};
%     end
%     contributions = cat(1, contributions{:});
%     feature_params.contributions = contributions;
    
    %% Feature extraction.
    [features]=featureExtractionDemo(exportArr, categoryArrIdx, poseArr, feature_params, []);
  
    %% Learn models.
    [models] = CategoryPoseLearningDemo( features, categoryArrIdx, poseArr);
%    [models] = CategoryPoseLearningDemo( features, categoryArrIdx, poseArr, feature_params.integration_levels);
    save([pwd '/Category_Pose/models/' datasetName '_models.mat'], 'models', 'feature_params');
end

