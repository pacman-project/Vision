%> Name: featureExtractionDemo
%>
%> Description: This script extracts features from part realizations that are
%> detected by LHOP.
%>
%>
%> @param exporArr exported array of realizations
%> @param categoryArrIdx Category labels of each image (as indices).
%> @param poseArr Pose labels of each image.
%> @param feature_params parameters of feature extraction algorithms.
%> @param testImageSize Test image's dimensions (x,y).
%> 
%> @retval features 1 by feature_params.learnedLevel cell that contain features that are extracted at each feature_params.learnedLevel=1,...,feature_params.learnedLevel.
%> @retval category_labels Category labels of samples.
%>
%>
%> Author: Mete
%>
%> Updates
%> Ver 1.0 on 13.04.2014
%> Ver 1.1 on 14.04.2014 (input/output changes)
function [features]=featureExtractionDemo(exportArr, categoryArrIdx, testImageSize)
    features = cell(size(categoryArrIdx,1),1);
    for sample_itr=1:size(categoryArrIdx,1)
        locations = exportArr(exportArr(:,5) == sample_itr,2:3);
        
        % Extract features.
        maskImg2=zeros(testImageSize);
        iddx = sub2ind(testImageSize, locations(:,1), locations(:,2));
        maskImg2(iddx)=1;
        feat = HOG(maskImg2)';
        features{sample_itr}=feat;
    end
    features = cat(1, features{:});
end