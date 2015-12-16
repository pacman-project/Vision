% [inds,cidx] = kmedioids(D,k)
%
% Performs k-mediods clustering; only requires a distance matrix D and
% number of clusters k.  Finds cluster assignments "inds" to minimize the
% following cost function: 

% sum(D(inds==i,inds==i),2), summed over i=1:k 
    
% Determining cluster assignments and cluster centers are both done in an
% efficient, vectorized way.  Cluster assignment is O(nk) and cluster
% centering is O(k*(max cluster size)^2)
%
% INPUTS
% D: nxn all-pairs distance matrix
% k: number of clusters
%
% OUTPUTS
% inds: nx1 vector of assignments of each sample to a cluster id
% cidx: kx1 vector of sample indices which make up the cluster centers

function [inds,cidx] = kmedioids(D, frequency, k)

    n = size(D,1);

    % randomly assign centers:
    cidx = randperm(n);
    cidx = sort(cidx(1:k));

    iter = 0;
    while 1
        inds = assign_pts_to_clusters(D,cidx);
        [cidx, energy_next] = update_centers(D,inds,k);

        if iter>0 && energy_next == energy
            break;
        end
        energy = energy_next;

        fprintf('iter: %04d, energy: %.02f\n',iter,energy)
        iter = iter+1;
    end

    function inds = assign_pts_to_clusters(D, cidx)
        S  = D(cidx, :);
        [~, inds] = min(S, [], 1);
    end

    function [cidx,energy] = update_centers(D, inds, k)
        energy = nan(k,1);
        cidx = zeros(k,1);
        for i = 1:k
           ind = find(inds==i);
           [energy(i),minind] = min(sum(D(ind,ind),2));
           cidx(i) = ind(minind);
        end
        energy = sum(energy);
    end
end    
