% this is to draw a view-independant part that is composed of several
% sub-parts

function is_ok = VisualizePart_VI(partsV, nCl, circleRad, vecLen)
    
    for i = 1:nCl
        
%         % compute normals
        q = partsV(i, 4:end);
        q = qnorm(q);
        normTemp = qvrot(q, [0,0,1]);
        normTemp = normTemp * vecLen;
        
%         normTemp = partsV(i, 4:6);
        
        quiver3(partsV(i, 1), partsV(i, 2), partsV(i, 3), normTemp(1), normTemp(2), normTemp(3), 'color', 'blue');
        hold on
        plotCircle3D([partsV(i, 1), partsV(i, 2), partsV(i, 3)],[normTemp(1), normTemp(2), normTemp(3)], circleRad(i));
        hold on
    end
    
%     xlim([-0.02 0.02])
%     ylim([-0.02 0.02])
%     zlim([-0.02 0.02])
    is_ok = true;
end

