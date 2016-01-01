% this is a trial function for arrayfun, bsxfun, pagefun

% suppose we have a set of points, P

function pagefun_arrayfunTrial()
%     numFrames = 10000;
%     numPoints = 70;
%     dim = 4;
%     
%     matrices = rand(dim,dim,numFrames, 'gpuArray'); % local frames of reference
%     pointsIn = rand(dim, numPoints, numFrames, 'gpuArray'); % 50 points for each frame of reference 
%     pointsOut = rand(dim, numPoints, numFrames, 'gpuArray'); % 50 points for each frame of reference
% 
%     %% baseline method
%     tic
%     for i = 1:numFrames
%         pointsOut(:, :, i) = matrices(:,:,i) * pointsIn(:,:,i);
%     end
%     toc
%     
%     %% pagefun-based method
%     tic
%         C = pagefun(@mtimes, matrices, pointsIn);
%         CC = gather(C);
%     toc
%     a = 2;
%%  % trial for a bsxfun

%     size = 1000;
%     a = rand(size, size, 'gpuArray');
%     b = rand(size, 1, 'gpuArray');
%     
%     tic
% %     for i = 1:10
%         d = bsxfun(@minus, a, b);
% %     end
%     toc;
%     
%     a = 2;

%% another trial for pagefun (inverse)
%     matr = rand(5,5, 1000, 'gpuArray');
%     
%     b = pagefun(@inv, matr);
%     bb = gather(b);
%     
%     a = 2;

%% this is one more trial (pairwise distances)
%     X = rand(3, 10^4, 'gpuArray');
%     Z = rand(3, 10^4, 'gpuArray'); 
%     
%     tic;
%     
%     XX = reshape(X,[1,3, 10000, 1]);
%     ZZ = reshape(Z,[1,3, 1, 10000]); 
% 
%     table = pagefun(@minus, XX, ZZ);
%     table = gather(table);
%     toc
%     
%     
%     a = 2;

% this is to try pagefun to convert sets of points to local frames of
% reference (number of points varies)

points = (:,:,1) = [1,2,3;2,3,4];

end






