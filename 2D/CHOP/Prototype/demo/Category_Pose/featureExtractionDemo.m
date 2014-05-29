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
function [features]=featureExtractionDemo(exportArr, categoryArrIdx, poseArr, feature_params, testImageSize)
    features = cell(size(poseArr,1),1);
    for sample_itr=1:size(poseArr,1)
        locations = exportArr(exportArr(:,5) == sample_itr,2:3);
        
        if ~feature_params.isTesting
            testImageSize = feature_params.imageSizes{categoryArrIdx(sample_itr)};
        end
        
        % Change here. We'll make all windows the same to extract features
        % correctly.
        if (~feature_params.isTesting && (categoryArrIdx(sample_itr)) < 3) || (feature_params.isTesting && ismember(180, testImageSize))
            testImageSize = feature_params.imageSizes{4};
            locations = locations + repmat([53, 70], size(locations,1), 1);
        end
        
        % Extract features.
        maskImg2=zeros(testImageSize);
        iddx = sub2ind(testImageSize, locations(:,1), locations(:,2));
        maskImg2(iddx)=1;
        feat = HOG(maskImg2)';
        features{sample_itr}=feat;
    end
    features = cat(1, features{:});
end



% function [features]=featureExtractionDemo(exportArr, categoryArrIdx, poseArr, feature_params, testImageSize)
%    
%     features = cell(1, feature_params.learnedLevel);
%     for level_id=1:feature_params.learnedLevel
%                 
%         level_idx= exportArr(:,4)==level_id ;
%         level_elements=exportArr(level_idx,:);
%                 
%         for sample_itr=1:size(poseArr,1)
%     
%             sample_elements = level_elements(:,5) == sample_itr;
%             locations = level_elements(sample_elements,2:3);
% 
%             % Extract features.
%             if feature_params.isTesting
%                 [feat] = featureExtraction( locations, feature_params.feature_type, testImageSize, feature_params );
%             else
%                 [feat] = featureExtraction( locations, feature_params.feature_type, feature_params.imageSizes{categoryArrIdx(sample_itr)}, feature_params );
%             end
% 
%             clear sample_elements
%             clear locations
% 
%             % Get features (hop and hog)
%             hops=feat.hop(:);
%             hogs=feat.hog(:);
%             h_hops=hist(hops,feature_params.histSize);
%             h_hogs=hist(hogs,feature_params.histSize);
% 
%             clear feat
%             f_temp=[h_hops h_hogs];
%             features{1,level_id}(sample_itr,:)=f_temp;
%             clear f_temp
%         
%         end
%         
%     end
% end