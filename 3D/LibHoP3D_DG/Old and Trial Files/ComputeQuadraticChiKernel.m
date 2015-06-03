function [distMatrix] = ComputeQuadraticChiKernel(PP, Q, A, m, is_GPU_USED)

    disp('Computing kernles ... ');
    
    if is_GPU_USED
        distMatrix = gpuArray.zeros(lenP, lenQ);
        P = gpuArray.zeros(lenQ, dim);
    else
        distMatrix = zeros(lenP, lenQ);
        P = zeros(lenQ, dim);
    end

    for i = 1:lenP

        P = repmat(PP(i, :), lenQ, 1);

        Z= (P+Q)*A;
        % 1 can be any number as Z_i==0 iff D_i=0
        Z(Z==0)= 1;
        Z= Z.^m;
        D = (P-Q)./Z;
        % max is redundant if A is positive-semidefinite
        distMatrix(i,:) = sqrt( max(diag(D*A*D'), 0));

        if mod(i, 100) == 0
            i
        end

    end
    
    if is_GPU_USED
        distMatrix = gather(distMatrix);
    end

end

