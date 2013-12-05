% matirix X (n*d) data samples
% matrix Z (nClusters*d) - cluster centroids
% distances - (n, nCluster) - matrix wuth distances

function [distances] = Isodata_distances(X, Z, n, nClusters, is_diagonal, want_disp)

    distances = zeros(n, nClusters); % initialization
    
    if want_disp
        str = ['computing distances from each point to each cluster ...'];
        disp(str);
    end    
    
    for ii = 1:nClusters
        curCentroid = Z(ii,:);
        b = repmat(curCentroid,n,1);
        diff = (X - b).^2;
        b = sqrt(sum(diff, 2));
        distances(:,ii) = b;
        if want_disp
            if mod(ii,50) == 0
                str = [num2str(ii), ' out of ', num2str(nClusters), ' clusters'];
                disp(str);
            end
        end
    end
    
    if (is_diagonal) % just make it low diagonal
        for i = 1: nClusters
            for j = i : nClusters
                distances(i,j) = 0;
            end
        end
    end
end
        