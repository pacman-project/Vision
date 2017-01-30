function [L C] = gmeans(data, maxk, alpha, check_func)

n = size(data, 1);
nC = 1;

if n < 2
     L = ones(n, 1);
     C = data;
     return;
end

C(1,:) = mean(data);
prevnC = 0;

while nC <= maxk
    %2
    incr_C = 0;
    [~, clusterStarts] = min(pdist2(C, data), [], 1);
    [L, ~, model] = kmeans_fast(data', clusterStarts);
    L = L';
    C = model.means';
    nC = size(C,1);
    if nC <= prevnC
         break;
    end
    prevnC = nC;
    validC = ones(nC,1) > 0;
    for i = 1:nC,
        %4
        if nC + incr_C >= maxk
             continue;
        end
        
        subdata = data(L==i,:);
        if isempty(subdata)
             validC(i) = 0;
             continue;
        end
        [t new_C] = look_gaussian(subdata, alpha, C(i,:), check_func);
        if (t ~= 1)
            %5
            C(i, :) = new_C(1, :);
            C = [C; new_C(2, :)];
            incr_C = incr_C + 1;
        end;
    end;
    if (incr_C == 0), break; end;
    C1 = C(1:nC,:);
    C1 = C1(validC, :);
    C2 = C((nC+1):end,:);
    C = cat(1, C1, C2);
    %disp(incr_C);
    nC = size(C, 1);
end;

function [t CC val]=look_gaussian(data, alpha, c, check_func)

if (size(data,1) == 1)
    t = 1; CC = [data; data];
    return;
end;

% choose two new centers, step 2
[~, lambda, s] = lmsvd(data, 1, []);
m = s' * sqrt(2 * lambda / pi);
CC(1, :) = c + m; CC(2, :) = c - m;
[~, clusterStarts] = min(pdist2(CC, data), [], 1);

%3
[~, ~, model] = kmeans_fast(data', clusterStarts);
CC2=CC;
CC = model.means';
if ~isequal(size(CC), size(CC2))
     t = 1;
     CC = [];
     return;
end
v = CC(1, :) - CC(2, :);

%4
nd = size(data, 1);
X = data*v' / norm(v);
%X = data*v';
[t, val] = check_func(X, alpha);
