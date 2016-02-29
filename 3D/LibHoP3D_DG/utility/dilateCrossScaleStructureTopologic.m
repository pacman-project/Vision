function crossScaleStructureOut = dilateCrossScaleStructureTopologic(crossScaleStructure, F, V, receptiveFieldRadTemp) %(crossScaleStructure(:, sc+1), F, V, receptiveFieldRadTemp, lenF);

    FaceEmpty = crossScaleStructure == 0;
    VertEmpty = F(:, FaceEmpty);
    VertEmpty = VertEmpty(:);

    U = geodesicDistances(V, F, VertEmpty, 3000);
    
    a= 2;
    
%     figure;
%     clf;
%     options.face_vertex_color = max(0, (receptiveFieldRadTemp - U));
%     plot_mesh(V, F, options);
%     colormap jet(256);
%     axis equal;
    
    VertEmpty1 = U <= receptiveFieldRadTemp;
    
    crossScaleStructureOut = VertEmpty1(F);
    crossScaleStructureOut = sum(crossScaleStructureOut, 1);
    crossScaleStructureOut(crossScaleStructureOut>0) = 1;
    crossScaleStructureOut = ~crossScaleStructureOut';
    
%     crossScaleStructureOut1 = find(crossScaleStructureOut);
%     VisualizeTriangulation(F(:, crossScaleStructureOut1), V);
    
    a = 2;
end