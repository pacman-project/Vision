
    nCl = 700;
    numCombs = 100000;
    
    tic

    is_sparse = false;
    disp('creating a matrix');

    if ~is_sparse
        matr = zeros(nCl, nCl, nCl);
    else
        parfor i = 1:nCl
    %         a = sparse(zeros(nCl, nCl));
    %         matr{i} = a;

            matr{i} = sparse(zeros(nCl, nCl));
        end
    end
toc;

tic


for i = 1:3
    SparsityTrial1(nCl, numCombs, matr, is_sparse); 
end

toc;