function KernelTrial()

    PP = [1,1,1,1,1];
    QQ = [1,2,0,1,1;   1,1,1,1,1;    1,2,4,1,1];
    
    AA = [  0.1465    0.5386    0.1239    0.2085    0.9479
            0.1891    0.6952    0.4904    0.5650    0.0821
            0.0427    0.4991    0.8530    0.6403    0.1057
            0.6352    0.5358    0.8739    0.4170    0.1420
            0.2819    0.4452    0.2703    0.2060    0.1665];
    
    dd = QC(PP, QQ, AA, 0.5)
    
    
    % this is to try pdist2 with this distance measure
    
    rbfKernel = @(X,Y) exp(-sigma .* pdist2(X,Y,@QC).^2);
    
    
    
   

    function distFun = QC(P,Q,A,m)
        P = repmat(P, size(Q, 1), 1);
        Z= (P+Q)*A;
        % 1 can be any number as Z_i==0 iff D_i=0
        Z(Z==0)= 1;
        Z= Z.^m;
        D= (P-Q)./Z;
        % max is redundant if A is positive-semidefinite
        distFun= diag(sqrt( max(D*A*D', 0)));
    end


end

