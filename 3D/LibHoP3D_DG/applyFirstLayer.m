% this is to apply the first layer to the images

function [] = applyFirstLayer(list_depth, list_mask, lenF, statistics1Layer, outFolder, sigma, sigmaKernelSize, ...
                                dxKernel, isErrosion, discRadius, is_guided, r_guided, eps,  is_mask_extended, maxExtThresh1, maxExtThresh2)

    % read the statistics 
    load(statistics1Layer); % 'cluster1Centres', 'cluster1Lengths', 'thresh', 'nClusters', 'dataSetNumber'
    
    isX = true;
    isY = true;
    isTrim = true;
    
    colormap lines;
    cmap = colormap;
    
    for i = 1:lenF
        
        I = imread(list_depth{i});
        I = I(:,:,1);
        
        imtool(I, [min(min(I)), max(max(I))]);
        
        if dataSetNumber == 2
            mask = imread(list_mask{i});
        else
            mask = [];
        end
        
        [I, Ix, Iy, mask, r, c, is_successfull] = preliminaryProcessing(I, mask, isErrosion, discRadius, isX, isY, ...
                    isTrim, dxKernel, sigmaKernelSize, sigma, is_guided, r_guided, eps,  is_mask_extended, maxExtThresh1, maxExtThresh2);
        
        % imtool(Ix, [-3, 3]);
        if ~is_successfull
            continue;
        end
        
        I1 = zeros(r,c,3);
        
        for j = 1:r
            fx = Ix(j,:);
            [nearestClusters, ~, ~, ~] = discretizeLine(fx, c, nClusters, cluster1Centres, cluster1Lengths, thresh);
            for jj = 1:c
                I1(j,jj, 1) = cmap(nearestClusters(jj),1);
                I1(j,jj, 2) = cmap(nearestClusters(jj),2);
                I1(j,jj, 3) = cmap(nearestClusters(jj),3);
            end
        end
        
        % imtool(I1);
        
        I2 = zeros(r,c,3);
        for j = 1:r
            fy = Iy(j,:);
            [nearestClusters, ~, ~, ~] = discretizeLine(fy, c, nClusters, cluster1Centres, cluster1Lengths, thresh);
            for jj = 1:c
                I2(j,jj, 1) = cmap(nearestClusters(jj),1);
                I2(j,jj, 2) = cmap(nearestClusters(jj),2);
                I2(j,jj, 3) = cmap(nearestClusters(jj),3);
            end
        end
        
        I1(:,:,1) = I1(:,:,1) .* mask;
        I1(:,:,2) = I1(:,:,2) .* mask;
        I1(:,:,3) = I1(:,:,3) .* mask;
        
        I2(:,:,1) = I2(:,:,1) .* mask;
        I2(:,:,2) = I2(:,:,2) .* mask;
        I2(:,:,3) = I2(:,:,3) .* mask;
        
        % imtool(I2);
        str1 = [outFolder, num2str(i), 'x.png'];
        str2 = [outFolder, num2str(i), 'y.png'];
        imwrite(I1, str1, 'png');
        imwrite(I2, str2, 'png');
                    
        % return;
        
        i
    end
    
    
end