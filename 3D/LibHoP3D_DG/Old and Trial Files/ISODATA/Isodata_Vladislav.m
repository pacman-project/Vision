% this is Vladislav's implementation of Isodata algorithm which is
% described here:
% https://github.com/himanshusingh/Machine-Learning-Tools-for-Opticks/wiki/ISODATA-Algorithm
% http://www.cs.umd.edu/~mount/Projects/ISODATA/ijcga07-isodata.pdf

% I also made some minor modifications to adapt algorithm to the input data
% and to the particular task

% parameters
% X - input data where 
% frequencies of each point
% rangeMin, rangeMax - range of the data (presumably 1:9)
% kInit - initial number of clusters
% nMin - minimal number of points to form the cluster
% iterations - maximal number of iterations
% sigmaMax - maximal standard deviation
% LMin minimum requested distance between two clusters
% pMax - maximum number of pairs to be mearged per iteration
% is_descrete - indicates if cluster centres are in discrete (integer), or
% continuous (float) space

%my parameter - percentUnclustered - percentage of the pixels which can be
%not clustered.

function [Z, sampleLables] = Isodata_Vladislav(X, frequencies, rangeMin, rangeMax ,kInit, nMin, iterations, sigmaMax, LMin, pMax, want_disp, is_descrete)
    
    sumF = sum(frequencies);

    range = 20;
    nClustersIter = zeros(iterations, 1);
    
    if want_disp
        disp('Isodata starts...');
    end

    [n, d] = size(X);
    % (1) randomly assign cluster centres
%     Z = rand(kInit, d); % centers of clusters
%     Z = Z * (rangeMax - rangeMin) + rangeMin; % Z*8 + 1
    
%   (1) Assign cluster centres in a much better way
    % take kInit points with maximal frequencies
    [~,IX] = sort(frequencies, 'descend');
    
    if kInit > n
        kInit = n;
    end
    IX = IX(1:kInit);
    Z = X(IX, :);
    
    nClusters = kInit;
    sampleLables = zeros(n,1);


    for i = 1:iterations
        
        
        if want_disp
            str = [num2str(i), ' iteration ...'];
            disp(str);
        end
        
%       (2) assign each point to the closest cluster center
        clusterFrequencies = zeros(nClusters, 1);
        clusterLables = 1:1:nClusters;
        distances = Isodata_distances(X, Z, n, nClusters, false, false);     % distances from each sample to each cluster center

        % find minimum distance and mark all samples correspondingly
        minDists = min(distances, [], 2);
        for ii = 1:n
            curLine = distances(ii,:);
            r = find(curLine == minDists(ii));
            ind = r(1);  % index of the closest cluster
            sampleLables(ii) = ind;
            clusterFrequencies(ind) = clusterFrequencies(ind) + frequencies(ii);
        end 
        clear('distances'); 

%       (3) remove cluster centres with frequency fewer than nMin  
        deletedClusterInds = find(clusterFrequencies<nMin);
        nDeleted = length(deletedClusterInds);
        
        sumF_deleted = sum(clusterFrequencies(deletedClusterInds));
        
        if (nDeleted ~= 0) % there is something to delete
            replaced_marks = zeros(1, nDeleted);
            [sampleLables, clusterLables, clusterFrequencies, Z, nClusters] = Isodata_deleteCluster(deletedClusterInds, nDeleted, replaced_marks, sampleLables, clusterLables, clusterFrequencies, Z, nClusters);                                                                                         
        end
        
        if want_disp && nDeleted > 0
            str = ['(3) ', num2str(nDeleted) ,' clusters removed. Remaining ', num2str(nClusters) ,' clusters' ];
            disp(str);
        end
        
%     (4) move each cluster centers to the centroid  of the associated
%     % cluster points

        centroids = zeros(nClusters, d);
        sumPixels = zeros(nClusters, 1);
        for ii = 1:n
            curCluster = sampleLables(ii);
            if curCluster ~= 0 % some pixels are not marked!
                centroids(curCluster, :) = centroids(curCluster, :) + X(ii,:)* frequencies(ii);
                sumPixels(curCluster) = sumPixels(curCluster) + frequencies(ii);
            end
        end

        Z = centroids./repmat(sumPixels, [1,d]); % new centers of clusters are now there          

        if is_descrete
            Z = round(Z);
        end
      
      
%     (5) compute average distance for each cluster and overall average distance
      %
        [delta, deltaOverall] = Isodata_distances_step5(X, sampleLables, frequencies, clusterFrequencies, Z, n, nClusters);
        
%         if want_disp
%             str = ['(5) deltaOverall = ', num2str(deltaOverall)];
%             disp(str);
%         end

%     (6)
        bNeedDivide = 0; 
        bNeedJoint = 0; 
        if (i == iterations) 
            LMin = 0; 
            bNeedJoint = 1; 
        elseif nClusters<=kInit/range 
            bNeedDivide = 1; 
        elseif nClusters>=range*kInit 
            bNeedJoint = 1; 
        elseif mod(i,2)==0 
            bNeedJoint = 1; 
        else 
            bNeedDivide = 1; 
        end 
        
        if (bNeedDivide~=0)
%     (7) for each cluster compute standard deviation in each direction

            V = zeros(nClusters, d);
            denominators = zeros(nClusters, 1);
            for ii = 1:n
                curCluster = sampleLables(ii);
                if curCluster ~= 0 % some pixels are not marked!
                    curFrequency = frequencies(ii);
                    V(curCluster, :) =  V(curCluster, :) + curFrequency * (X(ii,:) - Z(curCluster,:)).^2;
                    denominators(curCluster) = denominators(curCluster) + curFrequency;
                end
            end
            
            V = V./repmat(denominators, [1,d]);
            V = sqrt(V);
%     (8) split big clusters
            is_split = 0;
            nClustersNew = nClusters;
            for ii = 1:nClusters
            % find maximal standard deviation and compare it to sigmaMax
                if max(V(ii,:)) > sigmaMax
                    if ((delta(ii) > deltaOverall && clusterFrequencies(ii)> 2*(nMin + 1)) ||(nClusters <= kInit/2))
                        % separate this cluster
                        ind = find(V(ii,:) == max(V(ii,:))); % dimension
                        ind = ind(1);
                        Z = [Z; Z(ii,:)];
                        Z(ii, ind) = Z(ii, ind) + 0.5 * V(ii,ind);
                        nClustersNew = nClustersNew + 1;
                        Z(nClustersNew, ind) = Z(ii, ind) - 0.5 * V(ii,ind);
                        is_split = is_split + 1; 
                    end
                end
            end
            nClusters = nClustersNew;
            if want_disp && is_split > 0
                str = ['(8). We have split ', num2str(is_split), ' clusters. Now we have ', num2str(nClusters), ' clusters'];
                disp(str);
            end

        end
        
        if (bNeedJoint~=0)

    %     (9) Compute pairvise interclass distances
            interClassDists = Isodata_distances(Z, Z, nClusters, nClusters, true, false); % in fact matrix is diagonal!
            % make it diagonal

    %     (10) % select the subset with intercluster distance less than LMin
            values = interClassDists(interClassDists <= LMin);
            values = values(values > 0);
            values = sort(values, 'ascend');
            len = length(values);
            if ~isempty(values)
                row = [];
                col = [];
                ii = 1;
                while ii <= len
                    [rowCur, colCur] = find(interClassDists == values(ii));
                    row = [row; rowCur];
                    col = [col; colCur];
                    ii = ii + length(rowCur);
                end

                % pMax - maximum number of pairs to be mearged per iteration
                % mearge clusters 
                ii = length(row); 
                deletedClusterInds = [];
                ndeleted = 0; 
                replaced_marks = [];
                while (ii > 0 && ndeleted < pMax)
                    % always mearge clusters row(1) -> col(1)
                    deletedClusterInds = [deletedClusterInds, row(1)];
                    replaced_marks = [replaced_marks, col(1)];
                    ndeleted = ndeleted + 1;

                    inds = find(row == row(1));
                    inds1 = find(col == col(1));
                    inds = [inds; inds1];
                    row(inds) = [];
                    col(inds) = [];

                    ii = length(row);
                end
                % recompute cluster centres and frequencies
                for ii = 1:ndeleted
                    x1 = deletedClusterInds(ii);
                    x2 = replaced_marks(ii); % cluster indexes
                    n1 = clusterFrequencies(x1);
                    n2 = clusterFrequencies(x2); % weights of both clusters
                    Z1 = Z(x1,:);
                    Z2 = Z(x2,:); % claster centres
                    clusterFrequencies(x1) = clusterFrequencies(x1) + clusterFrequencies(x2);
                    Z(x1,:) = (n1 * Z1 + n2 * Z2)/(n1 + n2); % new cluster centre (Z(x2,:) is to be deleted)
                end
                [sampleLables, clusterLables, clusterFrequencies, Z, nClusters] = Isodata_deleteCluster(deletedClusterInds, ndeleted, replaced_marks, sampleLables, clusterLables, clusterFrequencies, Z, nClusters);   
                if want_disp
                    str = ['(10). We have mearged ', num2str(ndeleted), ' pairs. Now we have ', num2str(nClusters), ' clusters'];
                    disp(str);
                end    
            end
        end       
        
%       (11) show statistics and go to step 2
        str = ['iteration = ', num2str(i)];
        disp(str);
        str = ['Number of clusters = ', num2str(nClusters)];
        nClustersIter(i) = nClusters;
        disp(str);
%         str = ['(5) deltaOverall = ', num2str(deltaOverall)];
%         disp(str);
        str = ['Percent of unclustered pixels = ', num2str(sumF_deleted/sumF)];
        disp(str);
    end
    plot(nClustersIter);
    a = 2;
    
end




% 
%     % (2) assign each point to the closest cluster center
%     distances = zeros(n, nClusters);
%     
%     for ii = 1:nClusters
%         curCentroid = Z(ii,:);
%         b = repmat(curCentroid,n,1);
%         d = X - b;
%         d = d.^2;
%         b = sqrt(sum(d, 2));
%         distances(:,ii) = b;
%     end
%     distances = sqrt(distances); 
%     % find minimum distance
%     minDists = min(distances, [], 2);
%     % assign
%     for ii = 1:n
%         curLine = distances(ii,:);
%         r = find(curline == minDist(ii));
%         clusterLables(ii) = r(1);
%         clusterFrequencies(r(1)) = clusterFrequencies(r(1)) + 1;
%     end
%     
%     % (3) remove cluster centres with fewer than nMin
%     for ii = 1:nClusters
%         if clusterFrequencies(ii) < nMin  
%             clusterLables(clusterLables == ii)
%             clusterFrequencies(ii) = 0;
%             deletedClusters = (deletedClusters, ii);
%         end
%     end
%     
%     % (4) move each cluster centers to the centroid  of the associated
%     % cluster points
%     
%     for ii = 1:kInit
%         cur
%     end



















