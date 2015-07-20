% This function computes some points withing the triangle

function [Vadd, Nadd, numPoints] = compute_AdditionalPoints(V, F, N, VCent, NCent, numPoints, maxNumPoints)

    lenF = length(F);
    Vadd = zeros(maxNumPoints + 1, 3, lenF);
    Nadd = zeros(maxNumPoints + 1, 3, lenF);

    numPoints = numPoints + 1;
    ids = numPoints > 0;
    Vadd(1, :, ids) =  VCent;
    Nadd(1, :, ids) =  NCent;
    
    for j = 2:maxNumPoints+1
        
        ids = numPoints == j;
        [w1,w2,w3] = rand3Weights(j);
        
        for i = 1:j 
            VaddTemp = w1(i)* V(:,F(1,ids)) + w2(i)*V(:,F(2,ids)) + w3(i)*V(:,F(3,ids));
            NaddTemp = w1(i)* N(:,F(1,ids)) + w2(i)*N(:,F(2,ids)) + w3(i)*N(:,F(3,ids));
            Vadd(i, :, ids) =  VaddTemp;
            Nadd(i, :, ids) =  NaddTemp;
        end
    end
    
    
    Vadd = permute(Vadd, [3, 2, 1]);
    Nadd = permute(Nadd, [3, 2, 1]);
    
%     % add the central point to the end of the list
%     numPoints = numPoints + 1;
%     ind = 1:1:lenF;
%     for i = 1:3 
%         Vadd(ind, i, numPoints) =  VCent(i, ind);
%         Nadd(ind, i, numPoints) =  NCent(i, ind);
%     end
    
end



function [w1,w2,w3] = rand3Weights(numPoints)

    a = rand(1, numPoints);
    b = rand(1, numPoints);
    w1 = min(a, b);
    w2 = abs(a - b);
    w3 = 1 - max(a, b);

end

