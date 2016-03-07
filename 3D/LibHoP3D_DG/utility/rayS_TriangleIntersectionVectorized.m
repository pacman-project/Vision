% this is to check intersection of all rays with the scene


function [t_allPos, t_allNeg, numTpos, numTneg] = rayS_TriangleIntersectionVectorized (o_all, d_all, p0, p1, p2, isGPU)

    epsilon = 0.00001;
    numRays = size(o_all, 1);
    
    numT = zeros(1, size(o_all, 1));

    if isGPU
        p0 = gpuArray(p0);
        o_all = gpuArray(o_all);
        d_all = gpuArray(d_all);
        p1 = gpuArray(p1);
        p2 = gpuArray(p2);
    end
    
    P0 = p0;
    P1 = p1;
    P2 = p2;

    parfor i = 1:numRays
        
        p0 = P0;
        e1 = P1-P0;
        e2 = P2-P0;
        
        o = o_all(i, :);
        d = d_all(i, :);
        
%       q = bsxfun(@cross, d', e2')';
        q = e2;
        q(:,1) = d(2).*e2(:,3) - d(3).*e2(:,2);
        q(:,2) = d(3).*e2(:,1) - d(1).*e2(:,3);
        q(:,3) = d(1).*e2(:,2) - d(2).*e2(:,1);

    %     a  = e1 * q';  %dot(e1,q); % determinant of the matrix M
        a = sum(e1 .* q, 2);
        ids1 = ~(a>-epsilon & a<epsilon);  %  flag = 0; u= 0; v = 0; t = 0;

        a = a(ids1);
        p0 = p0(ids1, :);
        q = q(ids1, :);
        e1 = e1(ids1, :);
        e2 = e2(ids1, :);

        f = a.^(-1);
        s = bsxfun(@plus, o, -p0);   % o -p0 ;
        u = f .* sum(s .* q, 2);   % f*(s*q');

        ids2 = u >= 0;

        e1 = e1(ids2, :);
        e2 = e2(ids2, :);
        s = s(ids2, :);
        f = f(ids2);
        u = u(ids2);

%         r = bsxfun(@cross, s, e1);
        r = s;
        r(:,1) = s(:,2).*e1(:,3) - s(:,3).*e1(:,2);
        r(:,2) = s(:,3).*e1(:,1) - s(:,1).*e1(:,3);
        r(:,3) = s(:,1).*e1(:,2) - s(:,2).*e1(:,1);
        
        v = f.* (r(:, 1) * d(1) + r(:, 2) * d(2) + r(:, 3) * d(3));  %dot(d,r);
        ids3  =  ~(v<0.0 | (u+v)>1.0);

        f = f(ids3);
        e2 = e2(ids3, :);
        r = r(ids3, :);

        t = f .* sum(e2 .* r, 2); %dot(e2,r); % verified!
        t = uniquetol(t, epsilon);
        tpos = t(t> epsilon);
        tneg = t(t< -epsilon);
        t_allPos{i} = tpos;
        t_allNeg{i} = tneg;
        numTpos(i) = length(tpos);
        numTneg(i) = length(tneg);

%         if mod(i, 2000) == 0
%             str = [num2str(i), ' out of ' num2str(numRays)];
%             disp(str);
%         end
    end
    
    if isGPU
        t_all = gather(t_all);
    end
end

