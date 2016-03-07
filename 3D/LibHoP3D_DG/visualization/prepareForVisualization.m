
function [V, F, likelihoodsOut] = prepareForVisualization(radii, Vout, Nout, likelihoods)

    if size(Vout, 1) ~= 3
        Vout = Vout';
    end
    if size(Nout, 1) ~= 3
        Nout = Nout';
    end

    disp('Preparing data for visualization');
    lenParts = length(radii);
%     V = zeros(4 * lenParts, 3);
%     F = zeros(2 * lenParts, 3);
%     likelihoodsOut = zeros(2*lenParts, 1);
    
    Vtemp = zeros(12, lenParts);
    Ftemp = zeros(6, lenParts);
    likelihoodsOutTemp = zeros(2, lenParts);
    
     
    for i = 1:length(radii)

        if radii(i) == 0
            continue;
        end

        position = Vout(:, i);
        Norm = Nout(:, i);
        Norm = Norm/norm(Norm);

        if abs(dot(Norm, [1;0;0])) < 0.3
            xtemp = [1;0;0];
        elseif abs(dot(Norm, [0;1;0])) < 0.3
            xtemp = [0;1;0];
        else
            xtemp = [0;0;1];
        end

        Ytemp = cross(Norm, xtemp);
        Xtemp = cross(Norm, Ytemp);
        Ytemp = radii(i)*(Ytemp/norm(Ytemp));
        Xtemp = radii(i)*(Xtemp/norm(Xtemp));
        
        %% not parallelizable
        
        % establish 4 vertices and two riangles here
%         V((i-1)*4 + 1, :) = (position + Xtemp + Ytemp)';  % vertex
%         V((i-1)*4 + 2, :) = (position + Xtemp - Ytemp)';  % vertex  
%         V((i-1)*4 + 3, :) = (position - Xtemp - Ytemp)';  % vertex  
%         V((i-1)*4 + 4, :) = (position - Xtemp + Ytemp)';  % vertex  
% 
%         F((i-1)*2 + 1, :) = [(i-1)*4 + 1, (i-1)*4 + 2, (i-1)*4 + 3]; % faces
%         F((i-1)*2 + 2, :) = [(i-1)*4 + 1, (i-1)*4 + 3, (i-1)*4 + 4];
%         
%         likelihoodsOut((i-1)*2 + 1) = likelihoods(i);
%         likelihoodsOut((i-1)*2 + 2) = likelihoods(i);

        %% speed up
        Vtemp(:,i) = [(position + Xtemp + Ytemp); (position + Xtemp - Ytemp); (position - Xtemp - Ytemp); (position - Xtemp + Ytemp)];
        Ftemp(:,i) = [(i-1)*4 + 1; (i-1)*4 + 2; (i-1)*4 + 3; (i-1)*4 + 1; (i-1)*4 + 3; (i-1)*4 + 4];
        likelihoodsOutTemp(:, i) = [likelihoods(i); likelihoods(i)];
    end
    
    V = reshape(Vtemp, 3, 4 * lenParts)';
    F = reshape(Ftemp, 3, 2 * lenParts)';
    likelihoodsOut = reshape(likelihoodsOutTemp, 1, 2* lenParts)';
end

