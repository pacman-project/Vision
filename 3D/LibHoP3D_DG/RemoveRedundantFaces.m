function [F] = RemoveRedundantFaces(V, F)
    
    lenV = size(V, 2);
    multiplier = lenV +1;
    F = sort(F, 1);
    score = F(1, :) * multiplier^2 + F(2, :) * multiplier + F(2, :) + F(3, :);
    
    [r,c] = hist(score, unique(score));
    a = 2;
    
end

