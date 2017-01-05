function [L C] = gmeans(data, maxk, alpha, check_func)
%1
if (nargin < 2), maxk = size(data, 1) / 8; end;
if (nargin < 3), alpha = 0; end;
if (nargin < 4), check_func = @AnDartest; end;

if (ischar(check_func)),
    if (strcmp(check_func, 'gamma')),
        check_func = @GamKStest;
    else
        check_func = @AnDartest;
    end;
end;

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

function [t CC]=look_gaussian(data, alpha, c, check_func)
if (alpha == 0), alpha = 0.05; end;

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
t = check_func(X, alpha);

function [test] = AnDartest(x,alpha)
x = (x - mean(x)) / std(x);
n = length(x);
if n < 3,
    test = 1;
    return,
else
    x = x(:);
    x = sort(x);
    fx = normcdf(x,mean(x),std(x));
    i = 1:n;

    S = sum( (((2*i)-1)/n) *(log(fx)+log(1-fx(n+1-i))) );
    AD2 = -n-S;

    AD2a = AD2*(1 + 0.75/n + 2.25/(n^2)); %correction factor for small sample sizes: case normal

     if (AD2a >= 0.00 && AD2a < 0.200);
         P = 1 - exp(-13.436 + 101.14*AD2a - 223.73*AD2a^2);
     elseif (AD2a >= 0.200 && AD2a < 0.340);
         P = 1 - exp(-8.318 + 42.796*AD2a - 59.938*AD2a^2);
     elseif (AD2a >= 0.340 && AD2a < 0.600);
         P = exp(0.9177 - 4.279*AD2a - 1.38*AD2a^2);
     else
         P = exp(1.2937 - 5.709*AD2a + 0.0186*AD2a^2);
     end
end

if P < alpha % 0.01 alpha in these conditions.
   test=1;
else
   test=0;
end

return,

function [test] = GamKStest(x,alpha)
x=abs(x);
param = gamfit(x);
p2=gamcdf(x,param(1),param(2));
[H2,]=kstest(x,[x,p2],alpha);
test = ~H2;
return,