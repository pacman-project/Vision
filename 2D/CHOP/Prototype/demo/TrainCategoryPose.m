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
    feature_params.learnedLevel = 3;
    feature_params.integration_levels=[1,2,3];
    feature_params.isTesting = 0;
    feature_params.testImageSize = [];
    
    % Create an index-based category array.
    imageIds = exportArr(:,5);
    categoryArrIdx = zeros(numel(categoryArr), 1);
    categoryStrs = {'Mug', 'Bowl', 'Tool', 'Pan'};
    for categoryItr = 1:4
       samples = cellfun(@(x) strcmp(x, categoryStrs{categoryItr}), categoryArr);
       categoryArrIdx(samples) = categoryItr;
    end
    feature_params.imageSizes = {[180 220], [180 220], [285 360], [285 360]};
    
    %% Learn size distributions for every class in order to eliminate irrelevant windows in testing.
    minCounts = zeros(4,1);
    for categoryItr = 1:4
        samples = find(categoryArrIdx==categoryItr);
        numberOfSamples = numel(samples);
        realizationCounts = zeros(numberOfSamples,1);
        for sampleItr = 1:numberOfSamples
            realizationCounts(sampleItr) = nnz(imageIds == samples(sampleItr));
        end
        minCounts(categoryItr) = round(min(realizationCounts)*0.85);
    end
    feature_params.minCounts = minCounts;
    
    %% Feature extraction.
    [features]=featureExtractionDemo(exportArr, categoryArrIdx, poseArr, feature_params, []);
  
    %% Learn models.
    [models] = CategoryPoseLearningDemo( features, categoryArrIdx, poseArr, feature_params.integration_levels);
    save([pwd '/Category_Pose/models/' datasetName '_models.mat'], 'models', 'feature_params');
end

