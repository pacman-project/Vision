commonRoot = 'C:\Projects\Vladislav\'; 
root = [commonRoot, 'LibHoP3D_1/'];
addPaths(root);

[ depthPathDefault, depthPath, elPath] = getPathToData(2, commonRoot, root);

[list_depth, list_mask, ~, lenF] = extractFileListWashington(false, depthPath{1}, depthPathDefault, false, 1.0);

% load('Temp/depth_files.mat');
% lenF = length(list_depth);
isOk = zeros(1, lenF);

parfor i = 1:lenF
    
    try
        I = imread(list_depth{i});
        isOk(i) = 1;
    catch err
        isOk(i) = 0;
    end
    
end

a = 2;