function parForTrial
    dim = 1000;
    niter = 10^4
    map = zeros(dim, dim);
    f = @increaseMap;

    parfor i = 1:niter
        r = randi(dim);
        c = randi(dim);
        feval(f,r,c);     
    end
    
    a = 2;
    
    function increaseMap(r,c)
        map(r,c) = map(r,c) + 1;
    end
end