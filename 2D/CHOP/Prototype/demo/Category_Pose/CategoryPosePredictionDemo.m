%> Name: CategoryPoseLearning
%>
%> Description: This script extracts features from part realizations that are
%> detected by LHOP.
%>
%>
%> @param feature_set 1 by L cell of extracted features, where L is the 
%> number of layers.
%> 
%> @retval LearnedModels  Category and pose models learned.
%>
%> @param labels Used for measuring performance. Random numbers can be used
%> is labels are not available,
%>
%> @param integration_levels 1 by L vector whose variables are level ids  
%> which will be integrated.
%> 
%> @retval predictions predicted category and pose values.
%>
%> Author: Mete
%>
%> Updates
%> Ver 1.0 on 13.04.2014

function [predictions] = CategoryPosePredictionDemo( features, LearnedModels, labels)
%% Integrate features
% 
% features=[];
% for level_id=1:size(integration_levels,2)   
%     features=[features feature_set{1,integration_levels(level_id)}];   
% end
%% Prediction for classification
cmd = '-b 1';
%  [predicted_category_label, accuracy, prob_estimates] = svmpredict(labels_sub, features, LearnedModels.category_model, cmd);
[predicted_category_label,~, prob_estimates] = svmpredict(labels, features, LearnedModels.category_model,cmd);
 category_prediction.labels=predicted_category_label;
 category_prediction.prob_estimates = prob_estimates;
 
%% Prediction for pose estimation
 cmd = '-b 1';

 [predicted_pose_label, ~, prob_estimates] = svmpredict(rand(size(features,1),1), features, LearnedModels.pose_model{1,predicted_category_label},cmd);
 pose_prediction.labels=predicted_pose_label;
 pose_prediction.prob=prob_estimates;

 %% Predict pose interval


%   [pose_interval_predictions] = predict(LearnedModels.pose_model_interval{1,predicted_category_label},features);
    
%   predictions.pose_interval_predictions=pose_interval_predictions;
%   predictions.pose_interval_confidence=pose_interval_confidence;

 %%
 
 predictions.category_prediction=category_prediction;
 predictions.pose_prediction=pose_prediction;
end