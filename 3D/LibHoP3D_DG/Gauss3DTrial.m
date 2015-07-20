function Gauss3DTrial

    vectColors = eye(3);
    vecLen = 2;
    
    is_transformation = true;
    boxSize = 3;

%     mu = [10,10,10];
%     sigma = [1,0,0; 0,2,0; 0,0,1]; 
%     V = mvnrnd(mu,sigma,1000);
%      V = V';

    fileName = 'D:\Input Data\Meshes\Aim@Shape_Selected_2.00\D00017.obj';
    [V, F, ~] = meshRead(fileName);
    V = V';
    
    if size(V,2) ~= 3
        V = V';
    end

    % find principal components and make the largest dimension equal to 1

    meanV = mean(V, 1);
    V = V - repmat(meanV, [size(V, 1), 1]);
    
%     scatter3(V(:, 1), V(:, 2), V(:, 3));
%     xlabel('x')
%     ylabel('y')
%     axis equal
%     
    [W, EvalueMatrix] = eig(cov(V));
    Evalues = diag(EvalueMatrix);
    [Evalues, ids] = sort(Evalues, 'descend');
     W = W(:, ids);
 
    if is_transformation
        W(:,1) = W(:,1) / norm(W(:,1));
        W(:,2) = W(:,2) / norm(W(:,2));
        W(:,3) = W(:,3) / norm(W(:,3));

        % homogenious coordinates
        W(:, 4) = [0;0;0];
        W(4,:) = [0,0,0,1];
        W = inv(W);

        V(:,4) = zeros(size(V, 1),1);
        V = (W * V')';
    end

    scatter3(V(:, 1), V(:, 2), V(:, 3));
    xlabel('x')
    ylabel('y')
    axis equal

%     if ~is_transformation
%         hold on
%         plotFrame(W, vecLen, vectColors, [0,0,0]);
% 
%         hold off
%     end

    xs = V(:, 1);
    ys = V(:, 2);
    zs = V(:, 3);

    lenX = max(xs) - min(xs)
    lenY = max(ys) - min(ys)
    lenZ = max(zs) - min(zs)

end

function plotFrame(V, vecLen, vectColors, point)
    for ii = 1:3
        curVect = V(:, ii);
        curColor = vectColors(ii, :);
        XX = [point(1)- vecLen * curVect(1), point(1) + vecLen * curVect(1)];
        YY = [point(2)- vecLen * curVect(2), point(2) + vecLen * curVect(2)];
        ZZ = [point(3)- vecLen * curVect(3), point(3) + vecLen * curVect(3)];

        plot3(XX, YY, ZZ, 'Color', curColor);
        hold on
    end
end