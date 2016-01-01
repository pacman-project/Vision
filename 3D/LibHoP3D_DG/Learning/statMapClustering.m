% this is to do clustering of the statistical maps (by ISODATA)
% 1 - Isodata 2 - k-medioids, 3 - hierarchical clustering
% 

function [pairsLeft, pairsRight] = statMapClustering(nPrevClusters, statMapProperties, pairClusteringOptions, offsetConventional, layerID)

    number_of_classes_unknown = 1;
    nCl = 7;
    load('Temp/statMap.mat');  % statMapRight, statMapLeft
    
    pairsLeft = [];
    pairsRight = [];
    
    for kk = 1:2 % left and right stat map
        if kk == 1 
            statMap = statMapLeft;
        elseif kk == 2
            statMap = statMapRight;
        end
            
        for i = 1:nPrevClusters % central subpart
            for j = 1:nPrevClusters  % left or right subpart
                if isempty(statMap{i,j})
                    continue;
                end
                if max(statMap{i,j}(:)) < pairClusteringOptions.sieveThresh
                    continue;
                end
                tempMatr = statMap{i,j};
                [pairsTemp, nCl] = statisticsClustering(tempMatr, statMapProperties, number_of_classes_unknown, nCl, pairClusteringOptions, kk, layerID, offsetConventional, i, j);
                str = ['central: ', num2str(i), ', left () or right: ', num2str(j)];
                disp(str);
                a = 2;
                curPairs = [i*ones(nCl, 1), j*ones(nCl, 1), pairsTemp];
                if kk == 1
                    pairsLeft = [pairsLeft; curPairs];
                elseif kk == 2
                    pairsRight = [pairsRight; curPairs];
                end
            end      
        end
    end
    
    allPlots = findall(0, 'Type', 'figure');
    delete(allPlots);

end

function [parts, nCl] = statisticsClustering(tempMatr, statMapProperties, number_of_classes_unknown, nCl, pairClusteringOptions, kk, layerID, offsetConventional, cID, nID)

    is_visualization = 1;
    is_visualization2 = 0;
    selCluster = 5;
    circleRad = 0.005;
    vectLen = 0.01;
    
    if layerID <= 4 
        lim = 0.02;
    elseif layerID <= 6
        lim = 0.06;
    end

    penaltyThresh = pairClusteringOptions.penaltyThresh;
    nCl_max = pairClusteringOptions.nCl_max;
    alpha = pairClusteringOptions.alpha;
    clusteringMethod = pairClusteringOptions.clusteringMethod;
    sieveThresh = pairClusteringOptions.sieveThresh;


    tempMatr = tempMatr - sieveThresh;
    tempMatr(tempMatr<0) = 0;

    %% Find non-zero elements in the statistical maps 
    ids = find(tempMatr>0);
    [i1,i2,i3,i4,i5,i6] = ind2sub(size(tempMatr), ids);
    points = [i1,i2,i3,i4,i5,i6];
    frequencies = tempMatr(ids);
    
    % clustring should be performed in world's units
    points = statMapToWorld(points, layerID, kk, statMapProperties, offsetConventional);
    

    nCl_max = min(nCl_max, size(points, 1));

    if size(points, 1) == 1 % only one non zero entry in the statistical map
        parts = points;
        nCl = 1;
        c = 1;
    else   % parform clustering from the statistical map
        D = Mixed_Eucl_Quat_dist(points, points, alpha);
        if clusteringMethod == 1
    %         kInit = round(nnz(tempMatr)/5);
    %         nMin = 100;
    %         iterations = 100;
    %         sigmaMax = 0.01;
    %         LMin = 0.004;
    %         pMax = 250;
    %         want_disp = true;
    %         is_descrete = false;
    %         
    %          [Z, sampleLables] = Isodata_QuaternionDist(points, frequencies, kInit, nMin, iterations, sigmaMax, LMin, pMax, want_disp, is_descrete);

        elseif clusteringMethod == 2
        %         [inds,cidx] = kmedioids(D, frequencies,  k);
        %         parts = points(cidx, :);

        elseif clusteringMethod ==3

            L = tril(D);
            d = L(L>0)';
            Z = linkage(d, 'weighted'); %'average');

            % trying to find the optimal number of clusters
            % distance from points to cluster centres should be SMALL
            % distance between cluster centres should be large enough

            if number_of_classes_unknown
                nCl = 0;
                curPenalty = 0;
                penalty = zeros(nCl_max, 1);

                while (nCl < nCl_max) && curPenalty < 10^(-5)
                    nCl = nCl + 1;
                    parts = zeros(nCl, 6);
                    c = cluster(Z,'maxclust', nCl);

                    for i = 1:nCl
                        ids = c == i;
                        temp = points(ids, :);
                        parts(i, :) = sum(temp,1) / sum(ids); % centre of each clusters
                    end

                    D = Mixed_Eucl_Quat_dist(parts, parts, alpha);
                    D = D + eye(size(D));
                    % check how many clusters are close to one another
                    temp = D(D<penaltyThresh);
                    temp = penaltyThresh - temp;
                    curPenalty = sum(temp);
                    penalty(nCl) = curPenalty;
                end

                if curPenalty > 10^(-5)
                    nCl = nCl - 1;
                end
            end

            parts = zeros(nCl, 6);
            c = cluster(Z,'maxclust', nCl);
            for i = 1:nCl
                ids = c == i;
                temp = points(ids, :);
                parts(i, :) = sum(temp,1) / sum(ids); % centre of each clusters
            end
        end
        
        if is_visualization
            partsV = parts;
            partsV(:, 4:6) = parts(:, 4:6);
            figure;

            VisualizePart(layerID-1, cID, [0,0,0], [0,0,1], [1,0,0], [0,1,0], false, circleRad, vectLen);
            hold on

            lim = 0.06;
            xlim([-lim lim])
            ylim([-lim lim])
            zlim([-lim lim])

            for i = 1:nCl
                VisualizePart(layerID-1, nID, partsV(i, 1:3), partsV(i, 4:6), [], [], true, circleRad, vectLen);
                hold on
            end
            
            
            a = 2;
        end
        
        if is_visualization2
            nCl = max(c);
            partsV = parts;
            partsV(:, 4:6) = parts(:, 4:6)/10000;
            
            temp1V = temp1;
            temp1V(:, 4:6) = temp1V(:, 4:6)/20000;
            
            figure;
            quiver3(0,0,0, 0,0,0.01);
            hold on
            
            xlim([-lim lim])
            ylim([-lim lim])
            zlim([-lim lim])
            
            plotCircle3D([0,0,0],[0,0,0.01], 0.005);
            for i = selCluster
                quiver3(partsV(i, 1), partsV(i, 2), partsV(i, 3), partsV(i, 4), partsV(i, 5), partsV(i, 6), 'color', 'red', 'LineWidth', 5);
                hold on
                plotCircle3D([partsV(i, 1), partsV(i, 2), partsV(i, 3)],[partsV(i, 4), partsV(i, 5), partsV(i, 6)], 0.005);
                hold on
            end
            
            np = round(size(temp1V, 1)/10);
            ids = randperm(size(temp1V, 1),np);
                        
            for i = 1:np
                quiver3(temp1V(ids(i), 1), temp1V(ids(i), 2), temp1V(ids(i), 3), temp1V(ids(i), 4), temp1V(ids(i), 5), temp1V(ids(i), 6), 'color', 'blue', 'LineWidth', 1);
                hold on
            end
            a = 2;
        end
                
        % normalize the vectors of the parts before returning them
        for i = 1:nCl
            parts(i, 4:6) =  parts(i, 4:6) / norm( parts(i, 4:6)); % centre of each clusters
        end
        
    end
end

% conversion from the units of the statistical map to the world's units
% kk shows if this is a right or a left pair
function points = statMapToWorld(points, layerID, kk, statMapProperties, offsetConventional)

    xyzStep = statMapProperties.xyzStep;
    vectStep = statMapProperties.vectStep;
    [~, ~, ~, centreXY, centreZ, centreAngle] = computeStatMapSizes(statMapProperties, offsetConventional);

    if layerID == 3 || layerID == 5
        if kk == 1
            adderX = offsetConventional / xyzStep;
        elseif kk == 2
            adderX = - offsetConventional / xyzStep;
        end
        adderY = 0;
    elseif layerID == 4 || layerID == 6
        if kk == 1
            adderY = offsetConventional / xyzStep;
        elseif kk == 2
            adderY = - offsetConventional / xyzStep;
        end
        adderX = 0;
    end

    % back conversion of the units
    points(:, 1) = points(:, 1) + adderX - centreXY;
    points(:, 2) = points(:, 2) + adderY - centreXY;
    points(:, 3) = points(:, 3) - centreZ;
    points(:, 1:3) = points(:, 1:3) * xyzStep;
    points(:, 4:6) = (points(:, 4:6) - centreAngle) * vectStep;
end

function plotCircle3D(center,normal,radius)

    theta=0:0.01:2*pi;
    v=null(normal);
    points=repmat(center',1,size(theta,2))+radius*(v(:,1)*cos(theta)+v(:,2)*sin(theta));
    plot3(points(1,:),points(2,:),points(3,:),'r-');
    hold on

end


%           quiver3(0,0,0, 0,0,0.01);
%           hold on
%           plotCircle3D([0,0,0],[0,0,0.01], circleRad);


%                 quiver3(partsV(i, 1), partsV(i, 2), partsV(i, 3), partsV(i, 4), partsV(i, 5), partsV(i, 6), 'color', 'blue');
%                 hold on
%                 plotCircle3D([partsV(i, 1), partsV(i, 2), partsV(i, 3)],[partsV(i, 4), partsV(i, 5), partsV(i, 6)], circleRad);
%                 hold on













