% measures distances between two sets of surfaces: X and Z, which are represented
% by vectors of the same length n

function [distances] = Integral_distances(X, Z, n, m, is_diagonal, want_disp)

    distances = zeros(n, m); % initialization
     
    
    for ii = 1:m
        curSurface = Z(ii,:);
        b = repmat(curSurface,n,1);
        
        diff = abs(X - b);
        b = sqrt(sum(diff, 2));
        distances(:,ii) = b;
    end
    
    if (is_diagonal) % just make it low diagonal
        for i = 1: m
            for j = i : m
                distances(i,j) = 0;
            end
        end
    end
end
        