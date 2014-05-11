% this is to compute curvedness of the element
% el is an 18-dimensional vector

function [ curvedness ] = compute4Curvedness(el)

    indX = 1:2:18;
    indY = 2:2:18;
    Xs = el(indX);
    Ys = el(indY);
    
    varX = var(Xs);
    varY = var(Ys);
    
    curvedness = varX + varY;
end

