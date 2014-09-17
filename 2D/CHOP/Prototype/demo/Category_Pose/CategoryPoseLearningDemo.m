%> Name: CategoryPoseLearningDemo
%>
%> Description: This script extracts features from part realizations that are
%> detected by LHOP.
%>
%>
%> @param feature_set 1 by L cell of extracted features, where L is the 
%> number of layers.
%> 
%> @param category_labels N by 1 vector, where N is the number of images and 
%> the variables are category labels of the images.
%>
%> @param pose_labels N by 1 vector, where N is the number of images and 
%> the variables are pose labels of the images.
%>
%> @param integration_levels 1 by L vector whose variables are level ids  
%> which will be integrated.
%> 
%> @retval LearnedModels Learned category and pose models.
%>
%> Author: Mete
%>
%> Updates
%> Ver 1.0 on 13.04.2014
function [LearnedModels] = CategoryPoseLearningDemo( features, category_labels, pose_labels)
    %% Training for classification
    cmd='-s 0 -t 0 -b 1';
    category_model = svmtrain(category_labels, features, cmd);

    %% Training for pose estimation

    num_classes=max(category_labels);
    %%
    pose_labels = round(pose_labels / 30);
    pose_labels(pose_labels == 12) = 0;
    pose_labels = pose_labels + 1;
    % Class-wise training
    pose_model_class_wise = cell(1, num_classes);
    for class_idx=1:num_classes
        cmd='-s 0 -t 0 -b 1';

        subset_id=find(category_labels==class_idx);

        f_subset=features(subset_id,:);
        p_subset=pose_labels(subset_id,:);

        pose_model = svmtrain(p_subset, f_subset, cmd);

        pose_model_class_wise{1,class_idx}=pose_model;
    end

    %%

    LearnedModels.category_model=category_model;
    LearnedModels.pose_model=pose_model_class_wise;
end