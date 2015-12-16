function [is_ok] = VisualizeTriangulation(F, V, likelihood, cmap, FaceAlpha, FaceColor, EdgeColor)

    if nargin < 7
        EdgeColor = 'none';
    end
    if nargin < 6
        FaceColor = 'flat';
    end
    if nargin < 5
        FaceAlpha = 0.5;
    end
    if nargin < 4
        load('colormapGray.mat');
    end
    if nargin < 3
        likelihood = ones(size(F, 2), 1);
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

