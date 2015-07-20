% upsamples and rotates meshes a such that their largest dimension becomes
% equal to dowsample_rate

function upsampleMeshes(list_input, lenF, dowsample_rate)

    parfor i = 1:lenF 
        [V, F, ~] = meshRead(list_input{i});
        if size(V,2) ~= 3
            V = V';
        end
        V = PCA_model(V);
        
%         scatter3(V(:, 1), V(:, 2), V(:, 3));
%         xlabel('x')
%         ylabel('y')
%         axis equal
        xs = V(:, 1);
        scale = max(xs) - min(xs);
        V = V * dowsample_rate/scale;
        write_mesh(list_input{i}, V, F);
    end

end

function V = PCA_model(V)

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
    
end

