function result = TestTest()
    n = 10^3;
    nCl = 10^2;

    A = gpuArray(rand(n,3));
    B = gpuArray(rand(nCl,3));

    [ii jj] = meshgrid(1:n,1:nCl);
    M = numel(ii);
    ii = gpuArray(ii);
    jj = gpuArray(jj);

    res = arrayfun(@(n) distvect( A(ii(n),:), B(jj(n),:) ), 1:M);
    result = 1;
end

function d = distvect(a,b)
    d = sqrt(sum((a-b).^2));
end