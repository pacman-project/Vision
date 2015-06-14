% this function computes the structure tensor of the surface, given a set
% of normals

% flat surface - one orientation
% edge - two orientations
% corner - three orientations

function [T, V, D] = structureTensor3(Ns)
    
    len = size(Ns, 1);
    T = zeros(3,3);
    
    parfor i = 1:len
        Norm = Ns(i, :);
        if norm(Norm) ~= 0
            Norm = Norm/norm(Norm);
            T = T + (Norm' * Norm)/len;
        end 
    end
    
    [V,D] = eig(T);
end

