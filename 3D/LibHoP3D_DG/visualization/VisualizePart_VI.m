% this is to draw a view-independant part that is composed of several
% sub-parts

function is_ok = VisualizePart_VI(partsV, nCl, circleRad)
    
    for i = 1:nCl
        quiver3(partsV(i, 1), partsV(i, 2), partsV(i, 3), partsV(i, 4), partsV(i, 5), partsV(i, 6), 'color', 'blue');
        hold on
        plotCircle3D([partsV(i, 1), partsV(i, 2), partsV(i, 3)],[partsV(i, 4), partsV(i, 5), partsV(i, 6)], circleRad);
        hold on
    end
    
%     xlim([-0.02 0.02])
%     ylim([-0.02 0.02])
%     zlim([-0.02 0.02])
    is_ok = true;
end

