%> Name: getInitialNodes
%>
%> Description: When given an image, this function is used to get initial
%> simple nodes by applying Gabor filters over edge mask of the image. Each
%> response peak is considered as a simple feature instance.
%>
%> @param img Input image
%> @param filterCount (default 6) Number of rotated filters to apply on the
%> edge map.
%> @param currentFolder The workspace folder.
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 18.11.2013

function [ nodes ] = getNodes( img, filterCount, currentFolder)
    %% Parameters are set here for now, but they will move to the config file.
    gaborFilterThr = 0.7;   % Maximum response * gaborFilterThr will be accepted as a valid response. In [0,1].
    
    %% Process each band separately with the gabor filters and output the nodes.
    nodes = [];
    img = rgb2gray(img);
    doubleImg = img;
    
    edgeImg = double(edge(doubleImg(:,:), 'canny'));

    % Run each filter on the band and collect responses
    for filtItr = 0:(filterCount-1)
        
       % If the filter does not exist, create it.
       if ~exist([currentFolder '/filters/filt' num2str(filtItr+1) '.png'], 'file') 
           theta = (pi/filterCount) * filtItr;
           gaborFilt = gabor_fn(1, theta, 2, 0, 0.3);

           % Normalize the filter to [0,1]
           gaborFilt = (gaborFilt - min(min(gaborFilt)))/(max(max(gaborFilt))-min(min(gaborFilt)));
           imwrite(gaborFilt, [currentFolder '/filters/filt' num2str(filtItr+1) '.png']);
       else
           gaborFilt = double(imread([currentFolder '/filters/filt' num2str(filtItr+1) '.png']));
       end

       % get response
       responseImg = conv2(edgeImg, gaborFilt, 'same');
       responseImg = responseImg/max(max(responseImg)) > gaborFilterThr;

       %% Consider each connected component in binary response as existence of a filter.
       conn = bwconncomp(responseImg, 4);
       centers = regionprops(conn, 'Centroid');

       prevNodeCount = size(nodes,1);
       nodes = [nodes; cell(numel(centers), 2)];
       for centerItr = 1:numel(centers)
            nodes(prevNodeCount + centerItr, :) = {(filtItr+1), centers(centerItr).Centroid};
       end
    end
end

