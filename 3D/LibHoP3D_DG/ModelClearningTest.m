% This is for model clearning

% folderName = 'C:\Users\vvk201\Desktop\Pacman Vision Objects\Warehouse_Clean_Vladislav_full\Warehouse';

%'C:\Users\vvk201\Desktop\Pacman Vision Objects\Warehouse_Clean_Vladislav\Warehouse'; 
%'F:\Warehouse\Warehouse\Bottle\models';
folderName = 'D:\Input Data\Meshes\Aim@Shape_selected';

list_in = ReadImageNames(folderName);
lenF = length(list_in);
list_input = list_in; % list_input = {}; 
k = 0;

for i = 1:lenF
    str = list_in{i};
%     idsSP = strfind(str, ' ');
%     if ~isempty(idsSP)
%         str1 = str;
%         str1(idsSP) = '_';
%         disp(str1);
%         a = 2;
%         movefile(str, str1);
%     end

    
    if strfind(str, 'models') & strfind(str, '.obj')
        list_input{k+1} = str;
        k = k+1;
    end     
end

% %% this section is for remaning of the rendered folders
% 
% baseFolderOut = 'D:\PacmanVisionDataset\Views32\Dist_850\';
% baseFolderDepthIn =  'D:\RenderOutput2\depths\';
% baseFolderImagesIn = 'D:\RenderOutput2\images\'; 
% 
% lenF = k;
% categoryNamePrev = '';
% 
% for i = 1:lenF
%     curStr = list_input{i};
%     k = strfind(curStr, '\');
%     categoryName = curStr(k(end-2)+1:k(end - 1)-1);
%     
%     if ~strcmp(categoryName, categoryNamePrev)  % a new category starts
%         categoryNamePrev = categoryName;
%         categoryFolder = [baseFolderOut, categoryName, '\'];
%         if ~exist(categoryFolder, 'dir')
%             mkdir(categoryFolder);
%         end
%         objectID = 1;
%     end
%     folderDepthIn =  [baseFolderDepthIn,  'D_', num2str(i),'\'];
%     folderImagesIn = [baseFolderImagesIn, 'D_', num2str(i),'\'];
%     folderDepthOut = [categoryFolder, categoryName, '_', num2str(objectID), '/depth'];
%     folderImagesOut = [categoryFolder, categoryName, '_', num2str(objectID), '/images'];
%     
%     movefile(folderDepthIn,  folderDepthOut);
%     movefile(folderImagesIn, folderImagesOut);
%     
%     objectID = objectID + 1;
% end


% 
for i = 1:lenF
    
%     [V1,N1,F1,DF,PI] = simpleObjReader(list_input{i});
    disp(i);
    list_input{i} 
    [V,F,N] = read_obj(list_input{i});
    
%     scatter3(V(1,:), V(2,:), V(3,:));
%     axis equal;
%     
%     a = 2;
    
    %% translate to the origin 
%     V_centre = (max(V,[],2) + min(V,[], 2))/2; %sum(V, 2)./size(V, 2);
%     V = V - repmat(V_centre, [1, size(V, 2)]);
%     
     dirs = [1,2,3];

    %% make the largest dimension equal to 1
%        V = (PCA_model(V', dirs))';
%      V(3,:) = -V(3,:); % flip the z axis
    
%     scatter3(V(1,:), V(2,:), V(3,:));
%     axis equal;
    
    % put to the new centre
%     center = [max(V(1, :)) + min(V(1, :)); max(V(2, :)) + min(V(2, :)); max(V(3, :)) + min(V(3, :))]/2;
%     V = V - repmat(center, [1, size(V, 2)]);
%     
    scale = max(V(dirs(3),:) - min(V(dirs(3),:)));
    V = V/scale * 1; 
    V = V(dirs,:);

%     V=  V*0.6;

    figure;
    plot_mesh(V, F);
    figure;
    scatter3(V(1,:), V(2,:), V(3,:));
    axis equal;
    
%      pause(1.0);
% %     VisualizeTriangulation(F, V);

    if mod(i, 10) == 0
        allPlots = findall(0, 'Type', 'figure');
        delete(allPlots);
    end
    
% % 
%     str = list_input{i};
%     str = strrep(str, 'dae.obj', 'obj');
%     write_obj(str, V, F);

end











