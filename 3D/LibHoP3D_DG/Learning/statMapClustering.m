% this is to do clustering of the statistical maps (by ISODATA)
% 1 - Isodata 2 - k-medioids, 3 - hierarchical clustering
% 

function [pairsLeft, pairsRight] = statMapClustering(nPrevClusters, statMapProperties, pairClusteringOptions, offsetConventional, layerID, lenFiles)

    number_of_classes_unknown = 1;
    nCl = 7;
    load('Temp/statMap.mat');  % statMapRight, statMapLeft
    
    ReadVocabulary(layerID - 1);
    
    pairsLeft = [];
    pairsRight = [];
    
    [sizeXY, sizeZ, sizeAngle, ~,~,~] = computeStatMapSizes(statMapProperties, offsetConventional);
    sizeMatr = [sizeXY,sizeXY,sizeZ,sizeAngle,sizeAngle,sizeAngle,sizeAngle];
    
    pairClusteringOptions.sieveThresh = pairClusteringOptions.sieveThresh * log(lenFiles);
    
    if layerID == 8 || layerID == 10
        str = ['Temp/layerCurvatures', num2str(layerID - 1), '.mat'];
        load('Temp/layerCurvatures7.mat');
    end
    
    % the grammar rule must be the following:
    % at layers 8, 10, ect. we group together only those combinations
     % high curvature -> high curvature
     % low and moderate curvature -> low and moderate curvature
    
    for kk = 1:2 % left and right stat map
        if kk == 1 
            statMap = statMapLeft;
        elseif kk == 2
            statMap = statMapRight;
        end
            
        for i = 1:nPrevClusters % central subpart
            for j = 1:nPrevClusters  % left or right subpart
                
                if layerID == 8 || layerID == 10
                    % check the additional grammar rule
                    a = ismember(i, highCurvatures);
                    b = ismember(j, highCurvatures);
                    if a ~= b;
                        continue;
                    end
                end
                
                if isempty(statMap{i,j,1})
                    continue;
                end
                if max(statMap{i,j,2}) <= pairClusteringOptions.sieveThresh
                    continue;
                end
                [pairsTemp, nCl, is_ok] = statisticsClustering(statMap{i,j,1}, statMap{i,j,2}, sizeMatr, statMapProperties, number_of_classes_unknown, nCl, pairClusteringOptions, kk, layerID, offsetConventional, i, j);
                if nCl == 0 || ~is_ok
                    continue;
                end
                str = ['central: ', num2str(i), ', left () or right: ', num2str(j)];
                disp(str);
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

function [parts, nCl, is_ok] = statisticsClustering(idds, freq, sizeMatr, statMapProperties, number_of_classes_unknown, nCl, pairClusteringOptions, kk, layerID, offsetConventional, cID, nID)

    is_visualization = 0;
    is_visualization2 = 1;
    circleRad = 0.0033;
    vectLen = 0.01;
    
    if layerID <= 4 
        lim = 0.012;
    elseif layerID <= 6
        lim = 0.036;
    elseif layerID <= 8
        lim = 0.12;
    elseif layerID <= 10
        lim = 0.36;
    end
    is_ok = true;

    largeClusterThresh = pairClusteringOptions.penaltyThresh;
    nCl_max = pairClusteringOptions.nCl_max;
    alpha = pairClusteringOptions.alpha;
    clusteringMethod = pairClusteringOptions.clusteringMethod;
    sieveThresh = pairClusteringOptions.sieveThresh;
    
    xyzStep = statMapProperties.xyzStep;
    quaternionSteps = statMapProperties.quaternionSteps;
    
    penalty = 0.005;  % penalize large number of clusters in the statistical maps


    %% Find non-zero elements in the statistical maps 
%     ids = find(freq>sieveThresh1);
%     iddsA = idds(ids);        % above the threshold
%     frequencies = freq(ids);
%     
%     [i1,i2,i3,i4,i5,i6] = ind2sub(sizeMatr, iddsA);
%     points = [i1;i2;i3;i4;i5;i6]';
%     numPoints = size(points,1);
% 
%     % clustring should be performed in world's units
%     points = statMapToWorld(points, layerID, kk, statMapProperties, offsetConventional);
%     D = Mixed_Eucl_Quat_dist(points', points', alpha);
%     
%     % filtering of the statistical map
%     freqNew = zeros(size(frequencies));
%     smoothRad = penaltyThresh;
%     
%     for i = 1:numPoints
%         idsNeigh = find(D(i, :) < smoothRad & D(i, :)  > 0);
%         d = D(i, idsNeigh);
%         freqTemp = frequencies(idsNeigh);
%         weights = max(0, ones(1, length(idsNeigh)) - d/smoothRad);
%         freqNew(i) = frequencies(i) + sum(weights .* freqTemp);
%     end
%     freq(ids) = sqrt(freqNew);
    
    %% extract the points according to their updated frequencies
    freq = freq - sieveThresh;
    freq(freq<0) = 0;
    ids = freq>0;
    iddsA = idds(ids);
    freqSome = freq(ids);
    freqSome = freqSome + sieveThresh;
    [i1,i2,i3,i4,i5,i6,i7] = ind2sub(sizeMatr, iddsA);
    points = [i1;i2;i3;i4;i5;i6;i7]';
    points = statMapToWorld(points, layerID, kk, statMapProperties, offsetConventional);
    numPoints1 = size(points,1);
    
%   normalize all quaternions

    Norms = sqrt(points(:, 4).^2 + points(:, 5).^2 + points(:, 6).^2 + points(:, 7).^2);
    points(:, 4) = points(:, 4)./Norms;
    points(:, 5) = points(:, 5)./Norms;
    points(:, 6) = points(:, 6)./Norms;
    points(:, 7) = points(:, 7)./Norms;
    
    if numPoints1 == 0
        parts = [];
        nCl = 0;
        return
    end
    
    nCl_max = min(nCl_max, size(points, 1));

    if size(points, 1) == 1 % only one non zero entry in the statistical map
        parts = points;
        nCl = 1;
        c = 1;
        is_ok = false;
        
    else   % perform clustering from the statistical map
        D = Mixed_Eucl_Quat_dist(points', points', alpha);
        
        % set diagonal elements to 0
        iddds = sub2ind(size(D), 1:size(D,1), 1:size(D,1));
        D(iddds) = 0;
        
        if max(D(:)) <= largeClusterThresh % only one cluster is needed
            parts = points;
            nCl = 1;
            c = 1;
            number_of_classes_unknown = false;
        end
        
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

        elseif clusteringMethod == 3


            
            % extract ids of the lower triangular matrix
            tempM = ones(size(D));
            L = tril(tempM) - eye(size(D));
            idsM = find(L);
            d = D(idsM)';
            
%             L = tril(D);
%             d = L(L>0)';
            Z = linkage(d, 'weighted'); %'average');

            % trying to find the optimal number of clusters
            % distance from points to cluster centres should be SMALL
            % distance between cluster centres should be large enough

            if number_of_classes_unknown

                bouldinIndex = zeros(1, nCl_max);
                
                % check if only one cluster would work
                DdT = Mixed_Eucl_Quat_dist(points', points', alpha);
                if max(DdT(:)) < largeClusterThresh
                    nCl = 1;
                else
                    for nCl = 2:nCl_max
                        parts = zeros(nCl, 7);
                        c = cluster(Z, 'maxclust', nCl);
                        sigmas = zeros(1, nCl);  % average intra-cluster distances

                        % check if all clusters are small enough
                        isOk = true;
                        for i = 1:nCl
                            ids = c == i;
                            temp = points(ids, :);
                            % compute largest 
                            DdT = Mixed_Eucl_Quat_dist(temp', temp', alpha);

                            if max(DdT(:)) >= largeClusterThresh % at least one cluster is too big
                                isOk = false;
                                break;
                            end
                        end
                        if ~isOk
                            continue;
                        end

                        for i = 1:nCl
                            ids = c == i;
                            temp = points(ids, :);
                            freqSomeSome = freqSome(ids);
                            parts(i, :) = sum(temp,1) / sum(ids);                          % centre of each clusters
                            %   normalize the quaternion
                            parts(i, 4:7) = qnorm(parts(i, 4:7));
                            Dd = Mixed_Eucl_Quat_dist(parts(i, :)', temp', alpha);
                            sigmas(i) = sum(Dd .* freqSomeSome)/sum(freqSomeSome);   % sum(Dd)/size(temp,1);
                        end

                        Cs = Mixed_Eucl_Quat_dist(parts', parts', alpha); % inter-cluster distances
                        % between cluster centres or all points in each cluster???

                        sumDB = 0;
                        for i = 1:nCl
                            ratios = zeros(1, nCl);
                            for j = 1:nCl
                                if i == j
                                    continue;
                                end
                                ratios(j) = (sigmas(i) + sigmas(j))/(Cs(i,j));
                            end
                            sumDB = sumDB + max(ratios);
                        end
                        sumDB = sumDB/nCl;
                        bouldinIndex(nCl) = sumDB;
                    end

                    bouldinIndex(bouldinIndex == 0) = 100;
                    
                    % penalty for large numger of clusters
                    if layerID <= 6
                        temp = [0,0,0,0,0,     penalty:penalty:(nCl_max - 5)*penalty];
                    elseif layerID >= 7
                        temp = [0,0,0,0,0,0,0, penalty:penalty:(nCl_max - 7)*penalty];
                    end
                    temp = temp(1:length(bouldinIndex));                     
                    bouldinIndex = bouldinIndex + temp;


                    [~, nCl] = min(bouldinIndex);
                end

            end
            
   % parts should have the following format: muPos(3 numbers), sigmaPos(9
   % numbers), qOrient (4 numbers), SigmaOrient(16 numbers) ( 3+ 9 + 4 +16 = 32)

            parts = zeros(nCl, 32);
            c = cluster(Z,'maxclust', nCl);
            for i = 1:nCl
                ids = c == i;
                temp = points(ids, :);
                numPointsT = size(temp, 1);
                freqs = freqSome(ids);
                
                tempAll = zeros(sum(freqs), 7);
                cur = 1;
                for j = 1:numPointsT
                    tempAll(cur:cur+freqs(j)-1,:) = repmat(temp(j, :), [freqs(j), 1]);
                    cur = cur+freqs(j);
                end
                numPoints = size(tempAll, 1);
%                 tempAll(:, 1:3) = tempAll(:, 1:3) * 10000;
%                 temp(:, 1:3) = temp(:, 1:3) * 10000;
   
                if size(tempAll, 1) <= 16 
                    % we need at least 20 values
                    times = ceil(20 / size(tempAll, 1));
                    tempAll = repmat(tempAll, [times, 1]);
                end

                % this is done for regularization
                tempAll(:, 1:3) = tempAll(:, 1:3) + xyzStep * 0.7 *randn(size(tempAll(:, 1:3)));
                tempAll(:, 4:7) = tempAll(:, 4:7) + quaternionSteps * 0.7 * randn(size(tempAll(:, 4:7)));
                
                % compute Gaussians over position and orientation
                GMModel = fitgmdist(tempAll(:, 1:3),1, 'RegularizationValue', 10^-13);
                curMu = GMModel.mu;
                curSigma = GMModel.Sigma;
                
                % compute Gaussians over position and orientation
                GMModel = fitgmdist(tempAll(:, 4:7),1, 'RegularizationValue', 10^-13);
                curMuOrient = GMModel.mu;
                curSigmaOrient = GMModel.Sigma;
                
                % save the cluster parametrization
                parts(i, 1:3) = curMu; % centre of each clusters
                parts(i, 4:7) = curMuOrient;
                parts(i, 8:16) = curSigma(:);
                parts(i, 17:32) = curSigmaOrient(:);
                
                % test all cluster points
            end
        end
        
        % do check for all points

%         numPointsT = size(points, 1);
%         scoreGauss = zeros(3, numPointsT);
%         scoreFischer = zeros(3, numPointsT);
%    
%         for i = 1:nCl
%             curMu = parts(i, 1:3);
%             kappa = parts(i, 16);
%             mu = parts(i, 4:6);
%             curSigma = reshape(parts(i, 7:15), [3,3]);
%             invSigma = inv(curSigma);
%             vonMisesMultiplier = kappa/(2 * pi * (exp(kappa) - exp(-kappa)));
%             
%             for jj = 1:numPointsT
%                 scoreGauss(i,jj) = sqrt((points(jj, 1:3) - curMu) * invSigma * (points(jj, 1:3) - curMu)');
%                 scoreFischer(i,jj) = 1 - (vonMisesMultiplier * exp(kappa*points(jj, 4:6)*mu')/kappa) / ((vonMisesMultiplier * exp(kappa*mu*mu')/kappa));
%             end
%         end
%             
%         scoreAll = scoreGauss + 3 *scoreFischer;
%         [m,cEst] = min(scoreAll(1:3, :), [], 1);
%         sc = [scoreAll; scoreGauss; scoreFischer; c'; cEst; m];  % for testing

        
        if is_visualization
            partsV = parts;
            figure;

            VisualizePart(layerID-1, cID, [0,0,0], [], [0,0,1], [1,0,0], [0,1,0], false, circleRad, vectLen);
            hold on
            
            xlim([-lim lim])
            ylim([-lim lim])
            zlim([-lim lim])

            for i = 1:nCl
                VisualizePart(layerID-1, nID, partsV(i, 1:3), partsV(i, 4:7), [], [], [], true, circleRad, vectLen);
                hold on
            end
        end
        
%         if is_visualization2  % this is to visualize distribution in aech clusters
%             
%             for j = 1:nCl
%                 ids = cEst == j;
%                 partsV = points(ids, :);
%                 partsV(:, 4:6) =  partsV(:, 4:6)/100;
%                 
%                 figure;
%                 VisualizePart(layerID-1, cID, [0,0,0], [0,0,1], [1,0,0], [0,1,0], false, circleRad, vectLen);
%                 hold on
%                 colors = eye(3);
%                 
%                 selCluster = size(partsV, 1);
%                 for i = 1:selCluster
%                     quiver3(partsV(i, 1), partsV(i, 2), partsV(i, 3), partsV(i, 4), partsV(i, 5), partsV(i, 6), 'color', colors(j, :), 'LineWidth', 1);
%                     hold on
%                 end
%                 VisualizePart(layerID-1, nID, parts(j, 1:3), parts(j, 4:6), [], [], true, circleRad, vectLen);
%                 
%                 xlim([-lim lim])
%                 ylim([-lim lim])
%                 zlim([-lim lim])
%                 a = 2;
%             end
             
%             plotCircle3D([0,0,0],[0,0,0.01], 0.005);
%             for i = selCluster
%                 quiver3(partsV(i, 1), partsV(i, 2), partsV(i, 3), partsV(i, 4), partsV(i, 5), partsV(i, 6), 'color', 'red', 'LineWidth', 5);
%                 hold on
%                 plotCircle3D([partsV(i, 1), partsV(i, 2), partsV(i, 3)],[partsV(i, 4), partsV(i, 5), partsV(i, 6)], 0.005);
%                 hold on
%             end
%             
%             np = round(size(temp1V, 1)/10);
%             ids = randperm(size(temp1V, 1),np);
%                         
%             for i = 1:np
%                 quiver3(temp1V(ids(i), 1), temp1V(ids(i), 2), temp1V(ids(i), 3), temp1V(ids(i), 4), temp1V(ids(i), 5), temp1V(ids(i), 6), 'color', 'blue', 'LineWidth', 1);
%                 hold on
%             end
%             a = 2;
%         end
    end
    
%     % normalize the vectors of the parts before returning them
%     for i = 1:nCl
%         parts(i, 4:6) =  parts(i, 4:6) / norm( parts(i, 4:6)); % centre of each clusters
%     end  
end

% conversion from the units of the statistical map to the world's units
% kk shows if this is a right or a left pair
function points = statMapToWorld(points, layerID, kk, statMapProperties, offsetConventional)

    xyzStep = statMapProperties.xyzStep;
    quaternionSteps = statMapProperties.quaternionSteps;
    [~, ~, ~, centreXY, centreZ, centreAngle] = computeStatMapSizes(statMapProperties, offsetConventional);

    if layerID == 3 || layerID == 5 || layerID == 7 || layerID == 9
        if kk == 1
            adderX = offsetConventional / xyzStep;
        elseif kk == 2
            adderX = - offsetConventional / xyzStep;
        end
        adderY = 0;
    elseif layerID == 4 || layerID == 6 || layerID == 8 || layerID == 9
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
    points(:, 4:7) = (points(:, 4:7) - centreAngle) * quaternionSteps;
end

function plotCircle3D(center,normal,radius)

    theta=0:0.01:2*pi;
    v=null(normal);
    points=repmat(center',1,size(theta,2))+radius*(v(:,1)*cos(theta)+v(:,2)*sin(theta));
    plot3(points(1,:),points(2,:),points(3,:),'r-');
    hold on

end


%           Parameters of Fischer von Mises distribution
%                 sumX = sum(tempAll(:, 4:6), 1);
%                 muEst = sumX/norm(sumX);
%                 R = norm(sumX)/numPoints;
%                 
%                 deniminator = (1-R^2);
%                 if abs(deniminator) < 10^(-4)
%                     kappaEst = R*(3 - R^2)/deniminator;  % 3 means 3 dimensional space
%                 else
%                     kappaEst = 100;
%                 end
%                 if kappaEst > 100
%                     kappaEst = 100;
%                 end



%                 scoreGauss = zeros(1, numPointsT);
%                 scoreFisherVonMises = zeros(1, numPointsT);
%                 for jj = 1:numPointsT
%                     scoreGauss(jj) = gaussMultiplier * exp(-0.5 * (temp(jj, 1:3) - curMu) * invSigma * (temp(jj, 1:3) - curMu)');
%                 end


%           quiver3(0,0,0, 0,0,0.01);
%           hold on
%           plotCircle3D([0,0,0],[0,0,0.01], circleRad);


%                 quiver3(partsV(i, 1), partsV(i, 2), partsV(i, 3), partsV(i, 4), partsV(i, 5), partsV(i, 6), 'color', 'blue');
%                 hold on
%                 plotCircle3D([partsV(i, 1), partsV(i, 2), partsV(i, 3)],[partsV(i, 4), partsV(i, 5), partsV(i, 6)], circleRad);
%                 hold on





%                 for nCl = 1:nCl_max
%  
% %                     parts = zeros(nCl, 6);
% %                     c = cluster(Z, 'maxclust', nCl);
% %                     for i = 1:nCl
% %                         ids = c == i;
% %                         temp = points(ids, :);
% %                         parts(i, :) = sum(temp,1) / sum(ids); % centre of each clusters
% %                         Dd = Mixed_Eucl_Quat_dist(parts(i, :)', temp', alpha);
% %                         withinClustDist(i) = withinClustDist(i) + sum(Dd)/size(temp,1);
% %                     end
% %                     tempD = Mixed_Eucl_Quat_dist(parts', parts', alpha);
% %                     betweenClusterDist(i) = sum(tempD(:))/(2*nCl);
%                     
% %                     for i = 1:nCl
% %                         ids = c == i;
% %                         temp = points(ids, :);
% %                         parts(i, :) = sum(temp,1) / sum(ids); % centre of each clusters
% %                     end
%                 end
                
%                 figure
%                 x = 1:nCl;
%                 plot(x, withinClustDist);
%                 hold on
%                 plot(x, betweenClusterDist);
                
%                 nCl_b = find(withinClustDist - pairClusteringOptions.b2wRatio * betweenClusterDist < 0);
%                 nCl = nCl_b(1);
                
%                 while (nCl < nCl_max) && curPenalty < 10^(-5)
%                     nCl = nCl + 1;
%                     parts = zeros(nCl, 6);
%                     c = cluster(Z,'maxclust', nCl);
% 
%                     for i = 1:nCl
%                         ids = c == i;
%                         temp = points(ids, :);
%                         parts(i, :) = sum(temp,1) / sum(ids); % centre of each clusters
%                     end
% 
%                     D = Mixed_Eucl_Quat_dist(parts, parts, alpha);
%                     D = D + eye(size(D));
%                     % check how many clusters are close to one another
%                     temp = D(D<penaltyThresh);
%                     temp = penaltyThresh - temp;
%                     curPenalty = sum(temp);
%                     penalty(nCl) = curPenalty;
%                 end
% 
%                 if curPenalty > 10^(-5)
%                     nCl = nCl - 1;
%                 end













