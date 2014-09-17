function min_labels = mec( data, varargin )
%MEC Minimum Conditional Entropy Clustering
%   labels = mec(X) partitions points in N-by-P data matrix X into 2
%   clusters. This partition minimizes entropy of posterior probabilities.
%   Result is N-by-1 labels matrix. labels start from 1 and contiguous.
%
%   [ ... ] = mec(..., 'PARAM1',val1, 'PARAM2',val2, ...) specifies
%   optional parameter name/value pairs to control the iterative algorithm
%   used by mec.
%
%   Parameters:
%       'c' - Maximum number of clusters. Default is 2. Result can contain
%       smaller number of clusters than this parameter, if MEC achieves a
%       lower entropy.
%
%       'alpha' - alpha in alpha-entropy (default 1 for Shannon's entropy)
%
%       'kernel'    - Kernel of Parzen window
%           0 : Hypercube (default)
%           1 : Gaussian
%
%       'width'     - Bandwidth of Parzen window (default is 1.5)
%               The bandwidth of i-th deimension is set to w*sigma, where
%               sigma is the estimated standard deviation of i-th
%               dimension.
%
%       'kmeans_i'  - Run <I> times k-means for initialization (default 3)
%
%       'mec_maxiter' - Algoritm does not iterate any further and
%           terminates with best solution found after number following this
%           parameter. (Default 20)
%
%       'obj_thres' - (Default 0.0000001)
%
%   Sample run:
%   mu1 = [1 2];
%   sigma1 = [3 .2; .2 2];
%   mu2 = [-1 -1];
%   sigma2 = [2 0; 0 1];
%   data = [mvnrnd(mu1,sigma1,1000);mvnrnd(mu2,sigma2,1000)];
% 
%   % this is equivalent to labels = mec(data);
%   labels = mec(data, 'c', 2, 'alpha', 1, 'kernel', 0, ...
%               'width', 1.5, 'kmeans_i', 3, 'mec_maxiter', 20, 'obj_thres', 0.0000001);
%   
%   % plot the results:
%   figure, plot(data(labels==1,1), data(labels==1,2), '.', 'LineWidth', 3), hold on;
%   plot(data(labels==2,1), data(labels==2,2), '.', 'LineWidth', 3, 'color', 'red');
%   hold off;
%
%   Ref: Haifeng Li, Keshu Zhang, Tao Jiang, "Minimum Entropy Clustering
%   and Applications to Gene Expression Analysis". 142-151 2004
%   conf/csb/2004 CSB
%

%   Revision History:
%   Rev 2 : Bug fix related to making labels contiguous, kmeans does not
%   treat empty clusters as error anymore
%   Rev 1 : Performance and memory optimization, additional parameters
%   Rev 0 : Initial Release

% *********************************************************************** %
% Important notice by CHOP developer: Although I have downloaded this file
% from the web a while ago, I was unable to find the license/file on the
% Internet as of 07.09.2014 despite the best of my efforts. If the authors
% would like me to add the license file, I would be very happy to do so.
% *********************************************************************** %

p = inputParser;
p.addRequired('data', @isreal );
p.addParamValue('c', 2, @(x) isreal(x) && x >= 2);
p.addParamValue('alpha', 2, @(x) isreal(x) && x > 0);
p.addParamValue('kernel', 0, @(x) any(ismember(x, [0 1])));
p.addParamValue('width', 1.5, @(x) isreal(x) && x > 0);
p.addParamValue('kmeans_i', 3, @(x) isreal(x) && x >= 1);
p.addParamValue('mec_i', 1, @(x) isreal(x) && x >= 1);
p.addParamValue('mec_maxiter', 20, @(x) isreal(x) && x >= 1);
p.addParamValue('obj_thres', 0.0000001, @(x) isreal(x) && x > 0);
p.parse(data, varargin{:});

conv_gri = @(result,row) logical([result(1:row-1), 0, result(row:end)]); % used to add sample itself to the list of its neighbors

data = p.Results.data;
alpha = p.Results.alpha;
kernel = uint8(p.Results.kernel);
width = p.Results.width;
kmeans_iter = uint32(p.Results.kmeans_i);
c = uint32(p.Results.c);
mec_i = uint32(p.Results.mec_i);
mec_maxiter = uint32(p.Results.mec_maxiter);
obj_thres = p.Results.obj_thres;
clear p;

epsilon = 0.0000001; % used at entropy calculation, when a variable is smaller than this its entropy is treated as zero

hypercube = kernel == 0;
gaussian = kernel == 1;

col_mean = mean(data);
data_submean = data - repmat(col_mean,size(data,1),1);
col_sqsum_mean = sum(data_submean .^ 2) / (size(data,1) - 1);
col_sqrootsum = (col_sqsum_mean) .^ 0.5;

h = width * col_sqrootsum;

% init neighbor
if gaussian
    w_dist = @(a,b,c) transpose(sum(transpose(4 * (((repmat(a,size(b,1),1)-b) .^2) ./ repmat((c .^ 2), size(b,1),1))))) <= 1.0;
    neighbor = pdist(data, @(X,Y) w_dist(X,Y, h));
elseif hypercube
    w_dist = @(a,b,c) all((abs(repmat(a,size(b,1),1)-b) ./ repmat(c,size(b,1),1)) <= 0.5, 2)';
    neighbor = pdist(data, @(X,Y) w_dist(X,Y,h));
end
plim = length(neighbor);
psz = uint32(round(((sqrt(1+8*plim)) + 1) / 2)); % size of square form (psz-by-psz)

% init cluster
kmeans_labels = kmeans( data, c, 'Replicates', kmeans_iter, 'Start', 'sample', 'EmptyAction', 'drop');
unique_l = unique(kmeans_labels(:));
new_labels = zeros(size(kmeans_labels), 'uint32');
for i = 1:length(unique_l) % ensures that kmeans labels are contiguous
    new_labels(kmeans_labels == unique_l(i)) = i;
end
kmeans_labels = new_labels;
clear new_labels unique_l;

for m = 1:mec_i
    labels = kmeans_labels;
    
    % init posterior
    most_neighbor_cluster = zeros(size(data,1), 1, 'uint16');
    k_j = zeros(size(data,1), c, 'uint16');
    for i = 1:size(data,1)
        nbors = conv_gri(neighbor(getrowind(i,psz,plim)),i);
        k_j(i,:) = histc(labels(nbors), 1:c)';
        most_neighbor_cluster(i) = find(k_j(i,:) == max(k_j(i,:)), 1);
    end
    
    prev_obj = Inf;
    min_obj = Inf;
    min_labels = [];
%    fprintf(1, 'MEC run %d:\n', m);
    iter = 1;
    while iter <= mec_maxiter
        perm = randperm(size(data,1));
        for i = 1:size(data,1)
            curr_sample = perm(i);
            old_entropy = 0;
            new_entropy = 0;
            if most_neighbor_cluster(curr_sample) ~= labels(curr_sample)
                curr_neighbors = find(conv_gri(neighbor(getrowind(curr_sample,psz,plim)),curr_sample));
                for j = 1:length(curr_neighbors)
                    n = curr_neighbors(j);

                    nnn = double(sum(k_j(n,:)));
                    scn = double(k_j(n,labels(curr_sample)));
                    mnc = double(k_j(n,most_neighbor_cluster(curr_sample)));

                    if alpha == 1
                        scn_avg = scn / nnn;
                        mnc_avg = mnc / nnn;
                        if scn_avg > epsilon
                            old_entropy = old_entropy + scn_avg * log(scn_avg);
                        end
                        if mnc_avg > epsilon
                            old_entropy = old_entropy + mnc_avg * log(mnc_avg);
                        end
                        scn_avg = (scn-1) / nnn;
                        mnc_avg = (mnc+1) / nnn;
                        if scn_avg > epsilon
                            new_entropy = new_entropy + scn_avg * log(scn_avg);
                        end
                        if mnc_avg > epsilon
                            new_entropy = new_entropy + mnc_avg * log(mnc_avg);
                        end
                    else
                        old_entropy = old_entropy + (scn / nnn) .^ alpha;
                        old_entropy = old_entropy + (mnc / nnn) .^ alpha;
                        new_entropy = new_entropy + ((scn-1) / nnn) .^ alpha;
                        new_entropy = new_entropy + ((mnc+1) / nnn) .^ alpha;
                    end
                end
                if alpha < 1.0
                    old_entropy = -old_entropy;
                    new_entropy = -new_entropy;
                end
                if new_entropy > old_entropy % in fact this is -new_entropy > -old_entropy
                    for k = 1:length(curr_neighbors)
                        n = curr_neighbors(k);
                        k_j(n,labels(curr_sample)) = k_j(n,labels(curr_sample))-1;
                        k_j(n,most_neighbor_cluster(curr_sample)) = k_j(n,most_neighbor_cluster(curr_sample))+1;
                        most_neighbor_cluster(n) = find(k_j(n,:) == max(k_j(n,:)),1);
                    end
                    labels(curr_sample) = most_neighbor_cluster(curr_sample);
                end
            end
        end

        if alpha == 1.0
            obj = 0;
            d = double(k_j) ./ repmat(sum(k_j,2),1,c);
            for i = 1:length(d(:))
                if d(i) > epsilon
                    obj = obj - d(i) * log(d(i));
                end
            end
        else
            obj = size(data,1);
            d = double(k_j) ./ repmat(sum(k_j,2),1,c);
            for i = 1:length(d(:))
                if d(i) > epsilon
                    obj = obj - d(i) ^ alpha;
                end
            end
        end
        obj = obj / size(data,1);
        if alpha < 1.0
            obj = -obj;
        end
 %       fprintf(1, 'MEC Iteration %d: \t%g\n', iter, obj);
        improve = prev_obj - obj;
        if improve < obj_thres && iter >= 5 % at least 5 iterations enforced
            break;
        else
            prev_obj = obj;
        end
        iter = iter + 1;
    end

    if obj < min_obj
        unique_l = unique(labels(:));
        new_labels = zeros(size(labels));
        for i = 1:length(unique_l) % ensures that labels are contiguous
            new_labels(labels == unique_l(i)) = i;
        end
        min_obj = obj;
        min_labels = new_labels;
    end
end

%fprintf(1,'\nMinimum objective func value: %g\n', min_obj);

end

function result = getrowind(row,sz,lim)
    % this function is used to get indices correspoding to square form of
    % vector returned by pdist. It provides fast row access without calling
    % function squareform and hence with smaller amount of memory.
    result = zeros(1,sz-1,'uint32');
    res_i = 1;
    start = uint32(row-1);
    tmp = uint32(sz-2);
    first_iter = row - 1;
    
    while first_iter > 0
        result(res_i) = start;
        start = start + tmp;
        tmp = tmp - 1;
        res_i = res_i + 1;
        first_iter = first_iter - 1;
    end
    
    if res_i == 1
        result = uint32((1:sz-1));
    else
        result(res_i:sz-1) = (result(res_i-1) + tmp + 2: lim - (sz-row)*(sz-row-1) / 2);
    end
end
