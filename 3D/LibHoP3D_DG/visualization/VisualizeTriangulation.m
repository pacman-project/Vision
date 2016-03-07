function [is_ok] = VisualizeTriangulation(F, V, likelihood, cmap, FaceAlpha, FaceColor, EdgeColor)

    if size(F, 2) ~= 3
        F = F';
    end
    if size(V, 2) ~= 3
        V = V';
    end

    if nargin < 7
        EdgeColor = 'none';
    end
    if nargin < 6
        FaceColor = 'flat';
    end
    if nargin < 5
        FaceAlpha = 0.1;
    end
    if nargin < 4
        load('colormapGray.mat');
    end
    if nargin < 3
        likelihood = ones(size(F, 1), 1);
    end
    
    ids = round(likelihood*64);
    ids = ids+1;
    ids(ids > 64) = 64;
    color = cmap(ids, :);
    
    trisurf(F,V(:,1),V(:,2),V(:,3), 'FaceColor', FaceColor, 'FaceVertexCData', color, 'FaceAlpha', FaceAlpha, 'EdgeColor', EdgeColor);
    
%     light('Position',[-1.0,-1.0,100.0],'Style','infinite');
    axis equal;
    lighting phong;
    
    is_ok = true;
end

