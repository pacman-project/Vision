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
function [ models ] = TrainCategoryModel( datasetName, exportArr, categoryArrIdx )
    %% Set learning parameters.
    feature_params.binSize=32;
    feature_params.cellSize=32;
	feature_params.numOrient=32;
    feature_params.histSize=64;
    feature_params.feature_type='all';
    feature_params.integration_levels=[1,2,3,4,5];
    feature_params.learnedLevel = max(feature_params.integration_levels);
    feature_params.isTesting = 0;
    feature_params.testImageSize = [32 32];
    feature_params.categoryLevel = 4;
    feature_params.imageSize = [32 32];
    
    %% Feature extraction.
    [features]=featureExtractionDemo(exportArr, categoryArrIdx, feature_params.imageSize);
  
    %% Learn models.
    [models] = CategoryPoseLearningDemo( features, categoryArrIdx, []);
    
    %% Save models.
    save([pwd '/demo/Category_Pose/models/' datasetName '_models.mat'], 'models', 'feature_params');
end

