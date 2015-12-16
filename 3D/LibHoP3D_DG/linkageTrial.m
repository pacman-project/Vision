function linkageTrial

    m = rand(4,3);
    d1 = pdist(m)
    d2 = pdist2(m,m);
    U = tril(d2);
    d = U(U>0)'
    
    a = 2;
end

