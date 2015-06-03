function [distMatrix] = ComputeChi_Squared_Kernels(PP, Q, is_GPU_USED)

    disp('Computing kernles ... ');
    
    [lenP, dim] = size(PP);
    lenQ = size(Q, 1);
    
    if is_GPU_USED
        distMatrix = gpuArray.zeros(lenP, lenQ);
        P = gpuArray.zeros(lenQ, dim);
    else
        distMatrix = zeros(lenP, lenQ);
        P = zeros(lenQ, dim);
    end

    for i = 1:lenP

        P = repmat(PP(i, :), lenQ, 1);

        Z = (P - Q).^2;
        Denom = P+Q;
        Denom(Denom == 0) = 1;
        Z = (Z./Denom);
        distMatrix(i,:) = sum(Z,2)';
        
        if mod(i, 100) == 0
            i
        end

    end
    
    if is_GPU_USED
        distMatrix = gather(distMatrix);
    end
end

