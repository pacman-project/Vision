function V = PCA_model(V, dirs)

    [W, EvalueMatrix] = eig(cov(V));
    Evalues = diag(EvalueMatrix);
    [Evalues, ids] = sort(Evalues, 'descend');
     W = W(:, ids);
 
    W(:,1) = W(:,1) / norm(W(:,1));
    W(:,2) = W(:,2) / norm(W(:,2));
    W(:,3) = W(:,3) / norm(W(:,3));

    % homogenious coordinates
    W(:, 4) = [0;0;0];
    W(4,:) = [0,0,0,1];
    W = inv(W);

    V(:,4) = zeros(size(V, 1),1);
    V = (W * V')';
    
    V = V(:,1:3);
%     V = V(:, dirs);
    
end
