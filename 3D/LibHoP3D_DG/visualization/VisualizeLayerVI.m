% this is to visualize the view-invariant layers

function is_ok = VisualizeLayerVI(layerID)


    if layerID == 3
        dd = load('Temp/Layer3/partSelection3.mat'); % 'partsOut', 'coverageOut', 'pairsAll', 'nNClusters';
        lim = 0.02;
    elseif layerID == 4
        dd = load('Temp/Layer4/partSelection4.mat'); % 'partsOut', 'coverageOut', 'pairsAll', 'nNClusters';
        lim = 0.02;
    elseif layerID == 5
        dd = load('Temp/Layer5/partSelection5.mat'); % 'partsOut', 'coverageOut', 'pairsAll', 'nNClusters';
        lim = 0.06;
    elseif layerID == 6
        dd = load('Temp/partSelection6.mat'); % 'partsOut', 'coverageOut', 'pairsAll', 'nNClusters';
        lim = 0.06;
    end
    lenOut = dd.nNClusters{layerID};

    vectLen = 0.01;
    circleRad = 0.005;
    computeAxis = false;
    Xtemp = [1,0,0];
    Ytemp = [0,1,0];
    Norm =  [0,0,1];
    position = [0,0,0];
    
    
    for j = 1:lenOut
        figure;
        is_ok = VisualizePart(layerID, j, position, Norm, Xtemp, Ytemp, computeAxis, circleRad, vectLen);
        xlim([-lim lim])
        ylim([-lim lim])
        zlim([-lim lim])
        a = 2;
    end

    is_ok = true;
end

%         partsV = [pairsAll(partsOut(j, 1), :); pairsAll(partsOut(j,3), :)];
%         partsV(:, 4:6) = partsV(:, 4:6)/100;
%         figure;
%         quiver3(0,0,0, 0,0,0.01);
%         hold on
% 
%         plotCircle3D([0,0,0],[0,0,0.01], 0.005);
%         hold on
% 
%         xlim([-0.02 0.02])
%         ylim([-0.02 0.02])
%         zlim([-0.02 0.02])
% 
%         for i = 1:nCl
%             quiver3(partsV(i, 1), partsV(i, 2), partsV(i, 3), partsV(i, 4), partsV(i, 5), partsV(i, 6), 'color', 'blue');
%             hold on
%             plotCircle3D([partsV(i, 1), partsV(i, 2), partsV(i, 3)],[partsV(i, 4), partsV(i, 5), partsV(i, 6)], 0.005);
%             hold on
%         end