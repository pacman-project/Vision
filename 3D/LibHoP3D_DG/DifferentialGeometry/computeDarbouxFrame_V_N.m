% estimate the local frame of reference at this point
% method of Gabriel Taubin (modified by Vladislav Kramarev)

function [V, values] = computeDarbouxFrame_V_N(Norm, V, N, curvThresh)  % dx, dy,

if nargin == 4
    curvThresh = 10^-3;
end

%% make sure input has correct dimension
if size(V, 1) > 3
    V = V';
end
if size(N, 1) > 3
    N = N';
end

xs = V(1,:);
ys = V(2,:);
zs = V(3,:);

if size(Norm,2) > 1
    Norm = Norm';
end
if xs(1) == 0 && ys(1) == 0 && zs(1) == 0
    xs = xs(2:end);
    ys = ys(2:end);
    zs = zs(2:end);
    N = N(:, 2:end);
end

Norm = Norm/norm(Norm);
zeroThresh = 10^-10;

%     % (1). compute normal of the point

%     tanX = [1, 0, dx];
%     tanY = [0, 1, dy];
%     Norm = cross(tanX, tanY)';
%     Norm = Norm/norm(Norm);


    % (2). estimate a matrix M
    numPoints = length(xs);
    
    % compute weights for each point (reverse to euclidian distances)
    dist = sqrt(xs.^2 + ys.^2 + zs.^2);  
    sd = sum(dist);
    w = dist/sd;

%% variation 2. Computation of curvatures based on normals (SLOW)
% 
%     for ii = 2:numPoints
% 
%         vj = [0,0,0]';
%         vi = V(:, ii);
%         ni = N(:, ii);
% 
%         vect = (vi - vj);
%         Tij = (eye(3) - Norm*Norm')*vect;
%         Tij = Tij / norm(Tij);
% 
%         % estimate directional curvature
%         kij = 2 * (1 - Norm' * ni) / (norm(vect)^2); 
%         MTemp = w(ii) * kij * Tij * Tij';
%         M = M + MTemp;
%     end
    
    %%    ATTEMP TO SPEED UP (variation 2)
    
    Tij = (eye(3) - Norm*Norm')*V;
    TijNorm = sqrt(sum(abs(Tij).^2,1));
    ViNormSquared = sum(abs(V).^2,1);
%     Tij = Tij ./ repmat(TijNorm, [3,1]);
    Tij = Tij * diag(sparse(1./TijNorm));
 
    %     angles = zeros(1, size(Tij, 2));
    leNN = size(Tij, 2);

    T1 = Tij(:,1);
    
    % compute angles (i.e. its sin and cos) for each Tij w.r.t Tij(:,1)
%     sinS = zeros(1, leNN);
%     for i = 1:length(sinS)
%         T2 = Tij(:, i);
%         cosS(i) = dot(T1, T2);
%         sinS(i) = det([T1,T2,Norm]);
%         angles(i) = atan2(sinT, cosT);
%     end
    

     % compute 3 minors of the matrix
    M1 = -det([T1(1), N(1); T1(2), N(2)]);
    M2 =  det([T1(1), N(1); T1(3), N(3)]);
    M3 = -det([T1(2), N(2); T1(3), N(3)]);
    
    Mins = [M3;M2;M1];
    
%     MinsRep = repmat(Mins, [1, leNN]);
%     sinS = sum(MinsRep.*Tij ,1);   
    sinS = sum(diag(sparse(Mins)) * Tij, 1);

%     T11 = repmat(T1,[1, leNN]);
%     cosS = sum(T11.*Tij,1);
    cosS = sum(diag(sparse(T1)) *Tij, 1); 
       
    kij = 2 * (1 - Norm' * N) ./ ViNormSquared;
%     Tijm = Tij .* repmat(w, [3,1]) .* repmat(kij, [3,1]); 
%     Tijm = Tij * diag(sparse(w)) * diag(sparse(kij));
%     M = Tijm * Tij';

    M = Tij * diag(sparse(w)) * diag(sparse(kij)) * Tij';


    %%    ATTEMP TO SPEED UP (variation 1)

%     vi = [xs; ys; zs];
%     Tij = (eye(3) - Norm*Norm')*vi;
%     TijNorm = sqrt(sum(abs(Tij).^2,1));
%     ViNormSquared = sum(abs(vi).^2,1);
%     Tij = Tij ./ repmat(TijNorm, [3,1]);
%     
% %     angles = zeros(1, size(Tij, 2));
%     sinS = zeros(1, size(Tij, 2));
%     cosS = zeros(1, size(Tij, 2));
%     T1 = Tij(:,1);
%     % compute the angle for each Tij w.r.t Tij(:,1)
%     for i = 1:length(sinS)
%         T2 = Tij(:, i);
%         cosS(i) = dot(T1, T2);
%         sinS(i) = det([T1,T2,Norm]);
% %         angles(i) = atan2(sinT, cosT);
%     end
% 
%     % estimate directional curvature
%     kij = 2 * Norm' * vi ./ ViNormSquared; 
%     Tijm = Tij .* repmat(w, [3,1]) .* repmat(kij, [3,1]);
%     M = Tijm * Tij';

    
 %%   

    [V,D] = eig(M);
    values = abs(diag(D));
    ids =  find(values>zeroThresh);
    
    if length(ids) >= 2

        idsN =  find(values<zeroThresh);
        valTemp = values;

        values(ids(1)) = 3*valTemp(ids(1)) - valTemp(ids(2));
        values(ids(2)) = 3*valTemp(ids(2)) - valTemp(ids(1));
        [~, idx] = sort(abs(values), 'ascend');

        temp = idx(2);
        idx(2) = idx(3);
        idx(3) = temp;
        values = values(idx);


        V1 = V(:, ids(1)); V2 = V(:,ids(2));
        TT = [V1, -V1, V2, -V2];
    %     anglesT = zeros(1,4);
        curvTAvg = zeros(1,4);
        aT = 0.12; % angle thresh

        for i = 1:4
            T2 = TT(:, i);
            cosT = T1'*T2;
            sinT = det([T1,T2,Norm]);
    %         anglesT(i) = atan2(sinT, cosT);
            % select subset of vector pointing at this direction
            idsA = find(sinS < sinT+aT & sinS > sinT-aT & cosS < cosT+aT & cosS > cosT-aT);
            if ~isempty(idsA)
                curvTAvg(i) = sum(kij(idsA))/length(idsA);
            end    
        end
        curvTAvg = abs(curvTAvg);
        maxDir = find(curvTAvg == max(curvTAvg));
        maxDir = maxDir(1);
        firstPrinciple = TT(:, maxDir);
        if maxDir < 3 % V1
            V = [V(:,idsN(1)), firstPrinciple, V2];
        else % V2
            V = [V(:,idsN(1)), firstPrinciple, V1];
        end
        values = [0, 3*valTemp(ids(2)) - valTemp(ids(1)), 3*valTemp(ids(1)) - valTemp(ids(2))];
        values(2:3) = sort(values(2:3), 'descend');
    end
    
%     if this is a planar patch
    if abs(values(1)) < curvThresh  && abs(values(2)) < curvThresh  && abs(values(3)) < curvThresh
        V(:,1) = Norm;
        Xtemp = [1,0,0]';
        V(:,3) = cross(V(:,1), Xtemp);
        V(:,2) = cross(V(:,3), V(:,1));
    else
        
% %         if the direction of normal is reverced
        if abs(sum(V(:,1)+ Norm)) < 10^-5
            V(:,1) = -V(:,1);
        end
        id2 = find(abs(V(:, 2)) == max(abs(V(:, 2))));
        if V(id2(1), 2) < 0
            V(:, 2) = -V(:, 2);
        end
        id3 = find(abs(V(:, 3)) == max(abs(V(:, 3))));
        if V(id3(1), 3) < 0
            V(:, 3) = -V(:, 3);
        end

% %         if this is the umbilic point
        if abs(abs(values(2)) - abs(values(3))) < curvThresh*0.7 
%             make principal directions aligned with the main axis
            Xtemp = [1,0,0]';
            V(:,3) = cross(V(:,1), Xtemp);
            V(:,2) = cross(V(:,3), V(:,1));
        end
    end
    
    V(:, 1) = V(:,1)/norm(V(:,1));
    V(:, 2) = V(:,2)/norm(V(:,2));
    V(:, 3) = V(:,3)/norm(V(:,3));
    

end

%     % compute directional curvature in all 4 directions 
%     V1 = V(:, ids(1));  V2 = V(:, ids(2));  V3 = -V1;  V4 = -V2;
%      
%     disTT1 = Tij - repmat(V1, [1, size(Tij,2)]);
%     disTT2 = Tij - repmat(V2, [1, size(Tij,2)]);
%     disTT3 = Tij - repmat(V3, [1, size(Tij,2)]);
%     disTT4 = Tij - repmat(V4, [1, size(Tij,2)]);
%     
%     disTTT1 = sum(abs(disTT1), 1);
%     disTTT2 = sum(abs(disTT2), 1);
%     disTTT3 = sum(abs(disTT3), 1);
%     disTTT4 = sum(abs(disTT4), 1);
%     
%     [~,sortIndex1] = sort(disTTT1, 'ascend');
%     [~,sortIndex2] = sort(disTTT2, 'ascend');
%     [~,sortIndex3] = sort(disTTT3, 'ascend');
%     [~,sortIndex4] = sort(disTTT4, 'ascend');
%     
%     sortIndex1 = sortIndex1(1: 1+round(length(sortIndex1)*0.1));
%     sortIndex2 = sortIndex2(1: 1+round(length(sortIndex2)*0.1));
%     sortIndex3 = sortIndex3(1: 1+round(length(sortIndex3)*0.1));
%     sortIndex4 = sortIndex4(1: 1+round(length(sortIndex4)*0.1));
%     
%     sum1 = 0; sum2 = 0; sum3 = 0; sum4 = 0;
%     if length(sortIndex1)
%         sum1 = sum(abs(kij(sortIndex1)))/length(sortIndex1);
%     end
%     if length(sortIndex2)
%         sum2 = sum(abs(kij(sortIndex2)))/length(sortIndex2);
%     end
%     if length(sortIndex3)
%         sum3 = sum(abs(kij(sortIndex1)))/length(sortIndex3);
%     end
%     if length(sortIndex4)
%         sum4 = sum(abs(kij(sortIndex2)))/length(sortIndex4);
%     end
% 
%     ms = max([sum1, sum2, sum3, sum4]);
%     if ms == sum1
%         V = [V(:, idsN(1)), V1, V2];
%     elseif ms == sum2
%         V = [V(:, idsN(1)), V2, V1];
%     elseif ms == sum3
%         V = [V(:, idsN(1)), V3, V2];
%     elseif ms == sum4
%         V = [V(:, idsN(1)), V4, V1];
%     end


%     [V,D] = eig(M);
%     values = diag(D);
%     [valuesSort, idx] = sort(abs(values), 'ascend');
%     values = values(idx);
%     V = V(:, idx);
%     
%     k1 = 3*values(3) - values(2);
%     k2 = 3*values(2) - values(3);
%     values = [values(1), k1, k2];

