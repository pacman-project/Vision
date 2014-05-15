%%
clear
clc
    
    %% Feature Extraction

    export_location='data/export.mat';
    load(expor_location);
 
    feature_params.binSize=32;
    feature_params.cellSize=32;
	feature_params.numOrient=32;
    feature_params.histSize=64;
    feature_params.feature_type='all';
    
    Level=4;
    
    [features, category_label]=featureExtractionDemo(exportArr, categoryArr, poseArr, Level, feature_params);
    
    %% Split training and testing
    load(export_location)
    
    num_training=2951;
    
    random_index=randperm(size(category_label,2),size(category_label,2));
    tr_index=random_index(1:num_training);
    te_index=random_index(num_training+1:size(category_label,2));
    
    %% Training
    
    integration_levels=[1,2,3];
    
    for level_id=1:size(features,2)
    
    feature_set_tr{1,level_id}=features{1,level_id}(tr_index,:);
    
    end
    
    category_labels_tr=category_label(tr_index)';
    pose_labels_tr=poseArr(tr_index);    
    [LearnedModels] = CategoryPoseLearningDemo( feature_set_tr, category_labels_tr, pose_labels_tr, integration_levels);
    
    %% Testing

    for level_id=1:size(features,2)
    
    feature_set_test{1,level_id}=features{1,level_id}(te_index,:);
    
    end
    
    category_labels_te=category_label(te_index)';
    pose_labels_te=poseArr(te_index);   
    
    % Prediction is performed for each test sample:
    clear predictions
    [predictions] = CategoryPosePredictionDemo( feature_set_test, LearnedModels, category_labels_te, integration_levels);
        
    predicted_pose=predictions.pose_prediction.labels;
    
    calibrated_pose=mod(pose_labels_te,360);
    
    predicted_category=predictions.category_prediction.labels;
    
    

    
    
    
    