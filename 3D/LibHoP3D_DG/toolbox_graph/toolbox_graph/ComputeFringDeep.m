function [fringAll] = ComputeFringDeep(F, iter, lenF) % fringPrev, fring, lenF

    % compute distance from each face centre to each face centre
%     F1= F(1,:); F2 = F(2,:); F3 = F(3,:);
%     V1 = V(:,F1); V2 = V(:,F2); V3 = V(:,F3);
%     Fpos = (V1+V2+V3)/3;
%     Fdist = pdist2(Fpos', Fpos');

    gpuThresh  = 10^3;
    if lenF < gpuThresh % compute on CPU
        fring = compute_face_ring(F);
        adj1 =  sparse(flist2Adjacency(fring));
        adj1 = gpuArray(adj1);
        subtracter = eye(size(adj1));
        fringAll{1} = fring;
        adjCur = adj1;
        for j = 2:iter
            adjCur = adjCur * adj1 + adj1;
            adjCur(adjCur >1) = 1;
            adjCur = adjCur - subtracter;
        end
        fringAll{iter} = adjmatrix2list(adjCur);
        
    else  % compute on GPU
        fring = compute_face_ring(F);
        adj1 =  flist2Adjacency(fring);
        adj1 = gpuArray(sparse(adj1));
        subtracter = gpuArray.speye(size(adj1));
        fringAll{1} = fring;
        adjCur = adj1;
        for j = 2:iter
            adjCur = adjCur * adj1 + adj1;
            adjCur = ceil(mtimes(adjCur, 0.01));  % adjCur(adjCur >1) = 1; - does not work in this version of Matlab                
            adjCur = adjCur - subtracter;
        end
        
        % This is to do on GPU as much as possible
        [r,c] = find(adjCur);
        [r, ids] = sort(r, 'ascend');
        c = c(ids);
        r = r'; c = c';
        r = gather(r); c = gather(c);
    
%         idsD = sub2ind(size(Fdist), r, c);
%         distF = Fdist(idsD);
        
        fringAll{iter} = adjmatrix2listVladislav(r,c);

    end
end

%      adj1 = [0     1     0     0     0     0     1     0
%      1     0     1     0     0     0     0     0
%      0     1     0     1     0     1     0     0
%      0     0     1     0     1     0     0     0
%      0     0     0     1     0     1     0     1
%      0     0     1     0     1     0     1     0
%      1     0     0     0     0     1     0     0
%      0     0     0     0     1     0     0     0];


%     fringOut = {};
%     for jj = 1:lenF
%         cur = fring{jj};
%         temp = cur;
%         for jjj = 1:length(cur)
%             temp = union(temp, fringPrev{cur(jjj)});
%         end
%         fringOut{jj} = temp;
%     end