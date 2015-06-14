function  ArrayFuncDistTrial


    g = gpuDevice(1);
    
    tic
    
    m = 0.5;
    lenP = 900;
    lenQ = 10000;
    dim = 1000;
    
    
    
    A = gpuArray.rand(dim);
    PP = 2 * gpuArray.rand(lenP, dim);
    QQ = 2 * gpuArray.rand(lenQ, dim);
    distMatrix = gpuArray.zeros(lenP, lenQ);
    P = gpuArray.zeros(lenQ, dim);
    


%     A = rand(dim) + eye(dim);
%     PP = 2 * rand(lenP, dim);
%     QQ = 2 * rand(lenQ, dim);
%     distMatrix = zeros(lenP, lenQ);
% %     distMatrix1 = zeros(lenP, lenQ);
%     P = zeros(lenQ, dim);


%     for i = 1:lenP
%         
%         P = PP(i, :);
%         for j = 1:lenQ
% 
%             Q = QQ(j, :);
%             Z= (P+Q)*A;
%             1 can be any number as Z_i==0 iff D_i=0
%             Z(Z==0)= 1;
%             Z= Z.^m;
%             D= (P-Q)./Z;
%             max is redundant if A is positive-semidefinite
%             distMatrix(i,j) = sqrt( max(D*A*D', 0));
%             
%         end
%         
%         if mod(i, 10) == 0
%             i
%         end
%         
%     end 
    
    for i = 1:lenP
        
        P = repmat(PP(i, :), lenQ, 1);

            Q = QQ;
            Z= (P+Q)*A;
            % 1 can be any number as Z_i==0 iff D_i=0
            Z(Z==0)= 1;
            Z= Z.^m;
            D= (P-Q)./Z;
            % max is redundant if A is positive-semidefinite
            distMatrix(i,:) = sqrt( max(diag(D*A*D'), 0));

        if mod(i, 100) == 0
            i
        end
        
    end
    
%     max(max(distMatrix1 - distMatrix))
    

    toc;
    
    reset(g);
    
   
end

