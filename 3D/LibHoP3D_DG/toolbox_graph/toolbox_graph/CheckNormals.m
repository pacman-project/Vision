function N = CheckNormals(V, F, N, NF)

    if size(V, 2) ~= 3
        V = V';
    end
    if size(F, 2) ~= 3
        F = F';
    end
    if size(N, 2) ~= 3
        N = N';
    end
    if size(NF, 2) ~= 3
        NF = NF';
    end
    
    lenV = size(V, 1);
    % check it for each triangle
    
    for i = 1:lenF

        times = RaySceneIntersect(V(i, :), N(i, :), V, F);
        if mod(times, 2) == 1
            N(i, :) = -N(i, :);
        end
        if mod(i, 100) == 0
            disp(i);
        end    
    end

end



%     for i = 1:lenF
%         N0 = NF(i, :);
%         N1 = N(F(i, 1),:);
%         N2 = N(F(i, 2),:);
%         N3 = N(F(i, 3),:);
%         
%         a = N0 * N1';
%         b = N1 * N2';
%         c = N2 * N3';
%         
%         if a <= 0 || b<= 0 || c<= 0 % there is inconsistency. Check ray to scene intersections
%             times1 = RaySceneIntersect(V(F(i, 1), :), N1, V, F);
%             times2 = RaySceneIntersect(V(F(i, 2), :), N2, V, F);
%             times3 = RaySceneIntersect(V(F(i, 3), :), N3, V, F);
%             
%             if mod(times1, 2) == 1
%                 N(F(i, 1),:) = -N(F(i, 1),:);
%             end
%             if mod(times2, 2) == 1
%                 N(F(i, 2),:) = -N(F(i, 2),:);
%             end
%             if mod(times3, 2) == 1
%                 N(F(i, 3),:) = -N(F(i, 3),:);
%             end
%         end
%     end

