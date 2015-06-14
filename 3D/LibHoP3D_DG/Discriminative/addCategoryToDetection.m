load('statisticsSieved_3_2_7', '-mat');
load('statisticsAggregated_3_2_7','-mat');
load('depth_files.mat');
load('partsSelectionResults_3_2_7_a3.mat');

layerID = 3;
lenF = size(list_depth);
lenCombs = size(X, 1);
lenStat = size(statistics, 1);

FileCategoryID = zeros(lenF, 1); 
FileModelID = zeros(lenF, 1);

categoryNamePrev = '';
categoryID = 0;
modelNamePrev = '';
modelID = 0;

partIDs = zeros(lenStat, 1);
categoryIds = zeros(lenStat, 1);
modelIds = zeros(lenStat, 1);

for i = 1:lenF
    str = list_depth{i};
    inds = strfind(str, '/');
    if isempty(inds)
        inds = strfind(str, '\');
    end
    i3 = inds(end-2);
    i2 = inds(end-1);
    i1 = inds(end);
    categoryName = str(i3+1:i2-1);
    modelName = str(i2+1:i1-1);
    
    
    if ~strcmp(categoryName, categoryNamePrev)
        categoryNamePrev = categoryName;
        categoryID = categoryID + 1; 
        modelID = 0;
    end
    if ~strcmp(modelName, modelNamePrev)
        modelNamePrev = modelName;
        modelID = modelID + 1;
    end
    
    FileCategoryID(i) = categoryID;
    FileModelID(i) = modelID;
end
    
%--------------------------------------------------------------------------

nPrevClusters = nNClusters{layerID - 1} + 1; % because of empty cells

table = zeros(nPrevClusters, nPrevClusters, nPrevClusters);

% triple(left, center, right)
for ii = 1:lenCombs
    triple(X(ii, 1), X(ii, 2), X(ii, 3)) = ii;
end

parfor ii = 1:lenStat

    partIDs(ii) = triple(statistics(ii, 2), statistics(ii, 1), statistics(ii, 4));
    categoryIds(ii) = FileCategoryID(outputCoords(ii, 1));
    modelIds(ii) = FileModelID(outputCoords(ii, 1));
    if mod(ii, 10000) == 0
        ii
    end

end

 save('ready.mat', 'partIDs', 'categoryIds', 'modelIds');
    

   









