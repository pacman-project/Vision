%This is my trial of different implementation of pairwise distances

function d = pairiseDist_Trial(X, Z)

    n= size(X,1);
    nCl = size(Z,1);
    
%% Method 1  (11 sec)
%     tic
%     d = zeros(n, nCl);
%     % version 1
%     for i = 1:n
%         for j = 1:nCl
%             d(i,j) = sqrt(sum((X(i,:)-Z(j,:)).^2));
%         end
%         if mod(i,1000) == 0
%             disp(i);
%         end
%     end
%     toc;

%% method 2 (2.6 sec)
%     tic
%     d = zeros(n, nCl);
%     % version 1
%     parfor i = 1:n
%         for j = 1:nCl
%             d(i,j) = sqrt(sum((X(i,:)-Z(j,:)).^2));
%         end
%         if mod(i,1000) == 0
%             disp(i);
%         end
%     end
%     toc;

% %% method 3 (0.3 sec)
%     tic
%     d = zeros(n, nCl);
%     % version 1
%     for i = 1:n
%         temp = X(i,:);
%         d(i, :) = (Z(:,1) - temp(1)).^2 + (Z(:,2) - temp(2)).^2 + (Z(:,3) - temp(3)).^2;    
%         if mod(i,1000) == 0
%             disp(i);
%         end
%     end
%     d = sqrt(d);
%     toc;

% %% method 4 (0.09 sec)
    tic
    d = zeros(n, nCl);
    for i = 1:nCl
        temp = Z(i,:);
        d(:, i) = (X(:,1) - temp(1)).^2 + (X(:,2) - temp(2)).^2 + (X(:,3) - temp(3)).^2;    
    end
    d = sqrt(d);
    toc;

 %% method 5 (arrayFun) (0.09 sec)
    tic
    
    Xg = gpuArray(X);
    Zg = gpuArray(Z);

    
    dg = pagefun(@minus, Xg, Zg);
    dd = gather(dg);
    toc;
    a = 2;

end

