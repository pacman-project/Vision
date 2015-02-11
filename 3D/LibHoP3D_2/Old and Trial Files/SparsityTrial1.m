

function [] = SparsityTrial1(nCl, numCombs, matr, is_sparse) 

%     nCl = 700;
%     numCombs = 100000;
% 
%     is_sparse = false;
%     disp('creating a matrix');
% 
%     if ~is_sparse
%         matr = zeros(nCl, nCl, nCl);
%     else
%         parfor i = 1:nCl
%     %         a = sparse(zeros(nCl, nCl));
%     %         matr{i} = a;
% 
%             matr{i} = sparse(zeros(nCl, nCl));
%         end
%     end



    % our combinations to try
    X = rand(numCombs, 3);
    X = ceil(X * nCl);

    lenOutCoord = 3*numCombs;
    outCoords = rand(lenOutCoord, 3);
    outCoords = ceil(outCoords * nCl);
    values = zeros(1, lenOutCoord);



    disp('filling a matrix');


    % write combination to the 3D matrix
    for i = 1:numCombs
        if ~is_sparse
            matr(X(i,1), X(i,2), X(i,3)) = i;
        else
            matr{X(i,1)} (X(i,2), X(i,3)) = i;
        end

    %     if mod(i,2000) == 0
    %         i
    %     end
    end



    disp('indexing of matrix elements');


    % write combination to the 3D matrix
    for i = 1:lenOutCoord
        if ~is_sparse
            values(i) = matr(outCoords(i,1), outCoords(i,2), outCoords(i,3));
        else
            values(i) = matr{outCoords(i,1)}(outCoords(i,2), outCoords(i,3));
        end

    %     if mod(i,2000) == 0
    %         i
    %     end
    end



end

