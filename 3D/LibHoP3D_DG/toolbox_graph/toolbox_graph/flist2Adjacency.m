% list 

function A = flist2Adjacency(fring, is_sparse)
    
    if nargin == 1
        is_sparse = false;
    end
    
    lenF = length(fring);
    if is_sparse
        A = sparse(zeros(lenF, lenF));
    else
        A = zeros(lenF, lenF);
    end
    
    for i = 1:lenF
        A(i, fring{i}) = 1;
    end
end

