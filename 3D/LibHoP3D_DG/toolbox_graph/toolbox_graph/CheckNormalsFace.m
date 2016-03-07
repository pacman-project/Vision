% nPos are positions of the normals

function NF = CheckNormalsFace(V, F, nPos, NF)

    if size(V, 2) ~= 3
        V = V';
    end
    if size(F, 2) ~= 3
        F = F';
    end
    if size(nPos, 2) ~= 3
        nPos = nPos';
    end
    if size(NF, 2) ~= 3
        NF = NF';
        isTransposed = true;
    end
    
    isVisualization = true;
    
    lenF = size(NF, 1);
    % check it for each triangle
    
    isGPU = false;
     
    [t1,t2, numT1, numT2] = rayS_TriangleIntersectionVectorized (nPos,  NF, V(F(:, 1), :), V(F(:, 2), :), V(F(:, 3), :), isGPU);
    
    idsRev1 = numT1 > 0 & numT2 == 0;
    idsRev2 = numT1 == 1 & (numT2 == 2 | numT2 == 4 | numT2 == 6 );
    ids = idsRev1 | idsRev2;
    
    NF(ids, :) = -NF(ids, :);
    
    if isVisualization
        figure;
        VisualizeTriangulation(F, V);
        hold on

        ids = find(ids);
        if length(ids) > 50
            idsS = randperm(length(ids), 50);
            ids = ids(idsS);
        end

        for k = 1:length(ids)
            o = nPos(ids(k), :);
            d = NF(ids(k), :);
            quiver3(o(1), o(2), o(3), d(1), d(2), d(3));
            hold on
            a = 2;
        end
    end

    if isTransposed
        NF = NF';
    end
    

end



%       TT = [numT1;numT2]';
%       save('Temp/TT.mat', 'TT');
% %         aa = load('Temp/TT.mat');
% %         TT = aa.TT;

%     C = unique(TT, 'rows');
%     for j = 1:size(C, 1)
%         curComb = C(j, :)
%         [~, ind] = ismember(TT, curComb ,'rows');
% 
%         ids = find(ind);
%         figure;
%         VisualizeTriangulation(F, V);
%         hold on
%         
%         if length(ids) > 100
%             idsS = randperm(length(ids), 100);
%             ids = ids(idsS);
%         end
%         for k = 1:length(ids)
%             o = nPos(ids(k), :);
%             d = NF(ids(k), :);
%             quiver3(o(1), o(2), o(3), d(1), d(2), d(3));
%         end
% 
%         a = 2;
%     end


%     if is_visualize %&& mod(numIntersects, 2) == 1  % visualize wrong normals
%         figure
%         VisualizeTriangulation(F, V);
%         hold on
%         quiver3(o(1), o(2), o(3), d(1), d(2), d(3));
%         hold on
%         points = zeros(numIntersects, 3);
%         for i = 1: numIntersects
%             points(i, :) = o + d * t(i);
%         end
%         scatter3(points(:, 1), points(:,2), points(:, 3));
%         a = 2;
%     end


%         times = RaySceneIntersect(nPos(i, :), NF(i, :), V, F);


















%         if mod(times, 2) == 1
%             NF(i, :) = -NF(i, :);
%         end
%         if mod(i, 1000) == 0
%             str = [num2str(i), ' out of ', num2str(lenF)];
%             disp(str);
%         end 
        

