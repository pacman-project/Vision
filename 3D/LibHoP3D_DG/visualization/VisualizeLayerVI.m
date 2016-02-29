% this is to visualize the view-invariant layers

function is_ok = VisualizeLayerVI(layerID, startPoint)

    if nargin == 1
        startPoint = 1;
    end

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
        dd = load('Temp/Layer6/partSelection6.mat'); % 'partsOut', 'coverageOut', 'pairsAll', 'nNClusters';
        lim = 0.06;
    elseif layerID == 7
        dd = load('Temp/Layer7/partSelection7.mat'); % 'partsOut', 'coverageOut', 'pairsAll', 'nNClusters';
        lim = 0.12;
    elseif layerID == 8
        dd = load('Temp/Layer8/partSelection8.mat'); % 'partsOut', 'coverageOut', 'pairsAll', 'nNClusters';
        lim = 0.12;
    elseif layerID == 9
        dd = load('Temp/Layer9/partSelection9.mat'); % 'partsOut', 'coverageOut', 'pairsAll', 'nNClusters';
        lim = 0.6;
    elseif layerID == 10
        dd = load('Temp/Layer10/partSelection10.mat'); % 'partsOut', 'coverageOut', 'pairsAll', 'nNClusters';
        lim = 0.6;
    end
    lenOut = dd.nNClusters{layerID};
    
    strLayer = ['Temp/OR_node_layer_', num2str(layerID), '.mat' ];
    aa = load(strLayer);
    ORTable = aa.ORTable;
    
    
    vectLen = 0.01;
    circleRad = 0.0033;
    computeAxis = false;
    Xtemp = [1,0,0];
    Ytemp = [0,1,0];
    Norm =  [0,0,1];
    position = [0,0,0];
    QNorm = computeQuaternion(Norm, Norm);
    
   
    numInCluster = zeros(1, max(ORTable));
    
    for j = startPoint : length(ORTable)  % max(ORTable) % 
        if (layerID == 4 || layerID == 6) && j == 1
            figure();
            is_ok = VisualizePart(layerID, j, position, QNorm, Norm, Xtemp, Ytemp, computeAxis, circleRad, vectLen); % VisualizeTriangulation
        else
            if layerID == 4 || layerID == 6
                curCluster = ORTable(j-1);
            else
                curCluster = ORTable(j);
            end
            figure('Position',[20+40*(curCluster), 600 - 10 * numInCluster(curCluster),500,300]);
            is_ok = VisualizePart(layerID, j, position, QNorm, Norm, Xtemp, Ytemp, computeAxis, circleRad, vectLen);
            numInCluster(curCluster) = numInCluster(curCluster) + 1;
        end
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