% This function computes some points withing the triangle

function [Vadd, Nadd] = compute_AdditionalPoints(V, F, N, numPoints, maxNumPoints)

    lenF = length(F);
    Vadd = zeros(maxNumPoints, 3, lenF);
    Nadd = zeros(maxNumPoints, 3, lenF);

    for j = 1:maxNumPoints
        
        ids = numPoints == j;
        [w1,w2,w3] = rand3Weights(j);
        
        for i = 1:j 
            VaddTemp = w1(i)* V(:,F(1,ids)) + w2(i)*V(:,F(2,ids)) + w3(i)*V(:,F(3,ids));
            NaddTemp = w1(i)* N(:,F(1,ids)) + w2(i)*N(:,F(2,ids)) + w3(i)*N(:,F(3,ids));
            Vadd(i, :, ids) =  VaddTemp;
            Nadd(i, :, ids) =  NaddTemp;
        end
  
        j
    end
    
    Vadd = permute(Vadd, [3, 2, 1]);
    Nadd = permute(Nadd, [3, 2, 1]);
end



function [w1,w2,w3] = rand3Weights(numPoints)

    a = rand(1, numPoints);
    b = rand(1, numPoints);
    w1 = min(a, b);
    w2 = abs(a - b);
    w3 = 1 - max(a, b);

end

