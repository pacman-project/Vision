
% % this is to try different things with pairwise distances



% 
iter = 1000;
mSize = 3000; 

tic

X = rand(mSize, 3);
Y = rand(mSize, 3);


% D1 = sqrt(bsxfun(@plus, dot(X,X,2), dot(Y,Y,2)') - 2 * (X * Y'));
% 
% D = pdist2(X, Y);
% a = 2;

tic
parfor i = 1:iter 
    D = pdist2(X, Y);
end  % Elapsed time is 13.493408 seconds.
toc  

tic

X = gpuArray(X);
Y = gpuArray(Y);

for i = 1:iter
    D1 = sqrt(bsxfun(@plus, sum(X.^2, 2), sum(Y.^2, 2)') - 2 * (X * Y'));
end

D2 = gather(D1);
toc

a = 2;