% this is to compute clusters (boxes) of the third layer
% this function is an approximation of the function ComputeClusters3Layer
% The diffreence is that here we compute bounding boxes

clear all;

stats = load('statistics/statistics5D_3layer_19_5000.mat');
stats = stats.statistics;

gridX = 19;
gridY = 19;
gridZ = 71;
gridCenterX = ceil(gridX/2);
gridCenterY = ceil(gridY/2);

thresh1 = 0.000012; % shows if we take this cluster or not
thresh2 = 0.00001;

% those 14 = 2 + 3 + 9 + 1 + 1 components are el, eln, mu (3 elements), covariance
% matrix Sigma, score and index of the position (not restricted by thresh, to do a thresholding later)


cluster3stats = zeros(10000, 16); 
elNum = 0;

nClusters = 9;

%   9 2 6    
%   5 1 3    
%   8 4 7    
% displacements must be in this order {top, left, bottom right}
displacements = [-5, 0; 5, 0]; % [0, -5; 5, 0; 0, 5; -5,0]; % 5,-5; 5,5; -5,5; -5,-5];  % should work

Num = 2; %4; %8;
sumAll = sum(sum(sum(sum(sum(stats)))));
stats = stats/sumAll;

for clusterX = 1:9
    for clusterY = 1:9
        % compute the current central element
        el = compute2elementIndex(clusterX, clusterY, nClusters); % (clusterX - 1) * nClusters + clusterY;
        
        for clusterXn = 1:9 % this is an element we are looking around
            for clusterYn = 1:9
                eln = compute2elementIndex(clusterXn, clusterYn, nClusters); %(clusterXn - 1) * nClusters + clusterYn;
                               
                for i = 1:Num  % 4 for top-bottom-left-right; 8 for diagonal elements;
                    
                    % extract the slice corresponding to the displacement
                    curDisp = displacements(i,:);
                    sliceCenter = [gridCenterX + curDisp(1), gridCenterY + curDisp(2)];
                    slice1 = stats(el, eln, sliceCenter(1) - 2 : sliceCenter(1) + 2, sliceCenter(2) - 2 : sliceCenter(2) + 2, :);
                    slice1 = squeeze(slice1);
                    
                    sumSlice1 = sum(sum(sum(slice1)));
                    if sumSlice1 < thresh1  % for almost empty slices
                        continue;                
                    end
                                        
                    % otherwise we find a maximum point in this slice
                    [xm, ym, zm] = ind2sub(size(slice1), find(slice1 == max(max(max(slice1))))); % this must be a center of the cluster
                    xm = xm(1); ym = ym(1); zm = zm(1); % to avoid ambiguity
                                       
                    % now take another (smaller) slice of size slSize centered in this point
                    slSize = [5,5,11];
                    halfSliceSize = [floor(slSize(1)/2), floor(slSize(2)/2), floor(slSize(3)/2)];
                    slCenter = [ceil(slSize(1)/2), ceil(slSize(2)/2), ceil(slSize(3)/2)];
                    
                    % now we have to address the case if zm close to the
                    % margin
                    sliceCenter3D = [sliceCenter(1) + (xm - slCenter(1)), sliceCenter(2) + (ym - slCenter(2)), zm];    
                    
                    if zm > halfSliceSize(3) & zm < gridZ - halfSliceSize(3);                 
                        slice = stats(el, eln, sliceCenter3D(1)-halfSliceSize(1) : sliceCenter3D(1)+halfSliceSize(1), sliceCenter3D(2)-halfSliceSize(2) : sliceCenter3D(2)+halfSliceSize(2), zm - halfSliceSize(3): zm + halfSliceSize(3));                       
                    elseif zm < halfSliceSize(3)
                        margin = zm + halfSliceSize(3);
                        if mod(margin, 2) == 0
                            margin = margin + 1;
                        end
                        slice = stats(el, eln, sliceCenter3D(1)-halfSliceSize(1) : sliceCenter3D(1)+halfSliceSize(1), sliceCenter3D(2)-halfSliceSize(2) : sliceCenter3D(2)+halfSliceSize(2), 1 : margin);
                        sliceCenter3D(3) = 1 + floor((margin - 1)/2);
                    elseif zm >gridZ - halfSliceSize(3)
                        margin = zm - halfSliceSize(3);
                        if mod(margin, 2) == 0
                            margin = margin - 1;
                        end
                        slice = stats(el, eln, sliceCenter3D(1)-halfSliceSize(1) : sliceCenter3D(1)+halfSliceSize(1), sliceCenter3D(2)-halfSliceSize(2) : sliceCenter3D(2)+halfSliceSize(2), margin : gridZ);
                        sliceCenter3D(3) = margin + floor((gridZ-margin)/2);
                    end
                    
                    slice = squeeze(slice);
                    score = sum(sum(sum(slice)));  
                    
                    if (score) < thresh2
                        continue;
                    end
                        
                    %cluster inhibition :)
                    %stats(el, eln, sliceCenter3D(1)-slDim(1) : sliceCenter3D(1)+slDim(1), sliceCenter3D(2)-slDim(2) : sliceCenter3D(2)+slDim(2), zm - slDim(3): zm + slDim(3)) = stats(el, eln, sliceCenter3D(1)-slDim(1) : sliceCenter3D(1)+slDim(1), sliceCenter3D(2)-slDim(2) : sliceCenter3D(2)+slDim(2), zm - slDim(3): zm + slDim(3)) * 0.2;
                    
                    % estimate mu and sigma of the 3D Gaussian
                    [mu, Sigma] = estimate3DGaussian(slice);
                    % now transfer coordinates of mu from Slice coordinate
                    % system to world one
                    mu = mu + sliceCenter3D;                  
                    
                    elNum = elNum + 1;
                    cluster3stats(elNum, 1) = el;
                    cluster3stats(elNum, 2) = eln;
                    cluster3stats(elNum, 3:5) = mu;
                    cluster3stats(elNum, 6:14) = Sigma(:)';
                    cluster3stats(elNum, 15) = score; 
                    cluster3stats(elNum, 16) = i + 1;

                end
            end
        end
        
    end
end

cluster3stats = cluster3stats(1:elNum,:);

save('statistics/pairs3layer.mat','cluster3stats');