function U = geodesicDistances(vertex, faces, I, niter)

    if nargin == 3
        niter = 300; 
    end

    
    n = size(vertex,2);

    dotp = @(u,v)sum(u.*v,1);
    R = @(u)reshape(u, [1 1 length(u)]);
    Inv1 = @(M,d)[M(2,2,:)./d -M(1,2,:)./d; -M(2,1,:)./d M(1,1,:)./d];
    Inv  = @(M)Inv1(M, M(1,1,:).*M(2,2,:) - M(1,2,:).*M(2,1,:));
    Mult = @(M,u)[M(1,1,:).*u(1,1,:) + M(1,2,:).*u(2,1,:);  M(2,1,:).*u(1,1,:) + M(2,2,:).*u(2,1,:)];

    W = ones(n,1);

    U = zeros(n,1);
    ii = [faces(1,:) faces(2,:) faces(3,:) ];
    jj = [faces(2,:) faces(3,:) faces(1,:) ];
    kk = [faces(3,:) faces(1,:) faces(2,:) ];

    x  = vertex(:,ii);
    x1 = vertex(:,jj) - x;
    x2 = vertex(:,kk) - x;

    C = [R(dotp(x1,x1)) R(dotp(x1,x2)); R(dotp(x2,x1)) R(dotp(x2,x2))];
    S = Inv(C);
    a = sum(sum(S));
    d11 = sqrt(dotp(x1,x1));
    d11 = d11(:);
    d22 = sqrt(dotp(x2,x2));
    d22 = d22(:);

    iter = 1;
    err = 2;
    errs = zeros(niter, 1);
    
%     struct = {};
%     for i = 1:n
%         struct{i} = find(ii == i);
%     end
    struct2 = [];
    lengths = zeros(1, n);
    for i = 1:n
        temp = find(ii == i);
        struct2 = [struct2, temp];
        lengths(i) = length(temp);
    end
    cs = cumsum(lengths);
    starts = [1, cs(1:end-1) + 1];
    ends = cs;
    U1 = zeros(n, 1);

    while iter <= niter && err > 0.001
        uj = U(jj);
        uk = U(kk);
        u = [R(uj); R(uk)];
        w = R( W(ii) );

        b = dotp( sum(S,2), u );
        c = dotp( Mult(S,u), u ) - w.^2;
        delta = max( b.^2 - a.*c, 0);
        d = (b + sqrt(delta) )./a;
        alpha = Mult( S, u - repmat(d, 2, 1) );

        J = find( alpha(1,1,:)>0 | alpha(2,1,:)>0 );
        d1 = d11.*w(:) + uj(:);
        d2 = d22.*w(:) + uk(:);
        d = d(:);
        d(J) = min(d1(J), d2(J));

     %  U1 = accumarray(ii', d, [n 1], @min);  % too slow :(
     
        dd = d(struct2); 
        for i = 1:n
            U1(i) = min(dd(starts(i):ends(i)));  % nearly 2 times faster
        end

        U1(U1==0) = Inf;
        U1(I) = 0;
        err = sum(abs(U(:) - U1(:)));
        U = U1;
        iter = iter + 1;
        errs(iter) = err;
    end
end

        

            
%             A = arrayfun(@(x) min(d(struct{x})), 1:n, 'uniformoutput',true);
%            A = arrayfun(@(x) min(d(ii == x)), 1:n, 'uniformoutput',false);  % good stuff but too slow  

