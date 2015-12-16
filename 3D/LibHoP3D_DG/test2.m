folder = 'C:\Users\vvk201\Desktop\ObjectsWarehouseAligned';

list = load_filelist(folder);

% for i = 1:length(list)
filename = 'C:\Users\vvk201\Desktop\ObjectsWarehouseAligned\D_258.obj';
[V, F, N] = read_obj(filename);


VisualizeTriangulation(F, V);
    
%     strS = list{i};
%     strD = list{i};
%     strD = strD(1:end-4);
%     movefile(strS, strD);
%     
%     disp(i);
% end

