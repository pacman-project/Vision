%> Name: featureExtraction
%>
%> Description: This script extracts features from part realizations that are
%> detected by LHOP.
%>
%>
%> @param locations 1 by 2 vector, where the first and second columns 
%> represent the positions in x and y coordinates.
%> 
%> @retval features 1 by D vector, where the first and second columns 
%> represent the positions in x and y coordinates.
%>
%> Author: Mete
%>
%> Updates
%> Ver 1.0 on 10.04.2014

function [features] = featureExtraction( locations, feature_type, im_size, feature_params )
%% Histogram of Parts
if strcmp(feature_type,'hop') || strcmp(feature_type,'all')
    hop=mod(atan2(locations(:,1),locations(:,2)) .* 180./pi, 360   );
    features.hop=hop;
end
%% Histogram of Oriented Gaussians
if strcmp(feature_type,'hog') || strcmp(feature_type,'all')

binSize=feature_params.binSize;
cellSize=feature_params.cellSize;
numOrient=feature_params.numOrient;
%%

    im1=zeros(im_size(1),im_size(2));
    im2=im1(:);
  
    linearInd = sub2ind(size(im1), locations(:,1), locations(:,2));
    im2(linearInd)=1;
    im_temp = vec2mat(im2,size(im1,1));
    im_temp=im_temp';

%%

    hog = vl_hog(im2single(im_temp), cellSize, 'variant', 'dalaltriggs','numOrientations',numOrient) ;
    imhog = vl_hog('render', hog, 'variant', 'dalaltriggs','numOrientations',numOrient) ;

    features.hog=hog;
   %%
%     hog = vl_hog(im2single(im_temp), cellSize, 'variant', 'dalaltriggs','numOrientations',cellSize) ;
%     imhog2 = vl_hog('render', hog, 'variant', 'dalaltriggs','numOrientations',cellSize) ;
%     imshow(imhog2);
end

%% Von Neumann Entropy
if strcmp(feature_type,'von_neumann')  

    aa=(locations*locations');
    b=sqrt(sum((locations.*locations),2));
    b=b*b';

    c=aa./b;
    
    clear aa
    clear b
    
    ac=abs(acos(c));
    di=sum(ac,2);
    
    clear c
    
 %   -----------------------------------------
 %   Computed weighted Graph Laplacian
  
    l1=1-(diag(ac)./di);
           
    l2=-ac./(sqrt(repmat(di',[size(di,1) 1] )).*(repmat(di,[1, size(di,1)])));
    
    l2(1:length(l2)+1:numel(l2))=l1;
        
    clear l1;
    clear ac;
    clear di;
      
    v = eig(l2);
    von_neumann=-sum(v.*log(v));

    clear l2;
    clear v;
    features.von_neumann=von_neumann;

end


