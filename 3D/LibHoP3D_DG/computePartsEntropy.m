function [partEntropy] = computePartsEntropy(list_depth, statisticsLayerSieved, statisticsLayerAggregated, nPrevClusters, dataSetNumber)

    if dataSetNumber == 2
        numCategories = 51;
    end
        
    
    is_sparse = false;
    if nPrevClusters > 700
        is_sparse = true;
    end

    a = load(statisticsLayerSieved);
    outputCoords = a.outputCoords;
    statistics = a.statistics;
    
    clear('a');
    
    load(statisticsLayerAggregated);      % X, etc.
    
    lenF = size(list_depth, 1);
    lenCombs = size(X, 1);
    lenStat = size(statistics, 1);
    
    if dataSetNumber ~= 2   % TO FIX
        partEntropy = 3*ones(lenCombs, 1);
        return;
    end

    FileCategoryID = zeros(lenF, 1); 
    FileModelID = zeros(lenF, 1);

    categoryNamePrev = '';
    categoryID = 0;
    modelNamePrev = '';
    modelID = 0;

    partIDs = zeros(lenStat, 1);
    categoryIds = zeros(lenStat, 1);
    modelIds = zeros(lenStat, 1);
    
    indsCol = [2,1,4];
    statistics = statistics(:, indsCol);
    
%% define category for each file       

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
    
%% define 'partIDs', 'categoryIds', 'modelIds'

    
    % triples: left, center, right
    if is_sparse
        parfor i = 1:nPrevClusters
            triples{i} = sparse(zeros(nPrevClusters, nPrevClusters));
        end
    else  % if is_sparse
        triples = zeros(nPrevClusters, nPrevClusters, nPrevClusters);
    end
    
    
    % fill a table triples
    if ~is_sparse
        for i = 1:lenCombs
            triples(X(i, 1), X(i, 2), X(i, 3)) = i;
        end
    else
        for i = 1:lenCombs
            triples{X(i, 1)}(X(i, 2), X(i, 3)) = i;
        end
    end

    parfor ii = 1:lenStat
        
        if ~is_sparse
            partIDs(ii) = triples(statistics(ii, 1), statistics(ii, 2), statistics(ii, 3));
        else
            partIDs(ii) = triples{statistics(ii, 1)}(statistics(ii, 2), statistics(ii, 3));
        end
            
        categoryIds(ii) = FileCategoryID(outputCoords(ii, 1));
        modelIds(ii) = FileModelID(outputCoords(ii, 1));
        if mod(ii, 10000) == 0
            ii
        end
    end

     save('Temp/ready.mat', 'partIDs', 'categoryIds', 'modelIds');
     
%% build a table

    table = zeros(lenCombs, numCategories);

    for i = 1:lenStat

        if partIDs(i) > 0
            table(partIDs(i), categoryIds(i)) = table(partIDs(i), categoryIds(i)) + 1;
        end

        if mod(i, 100000) == 0
            i
        end
    end
    
%%  Normalize table and compute entropy from histograms

    denominator = sum(table, 1);
    multiplier = max(denominator);
    denom = repmat(denominator, [lenCombs, 1]);

    table = round(multiplier*(table ./ denom));  % adjast the histogram

    partEntropy = zeros(lenCombs,1);

    parfor i = 1:lenCombs
         partEntropy(i) = theirEntropy(table(i, :));

         if mod(i, 10) == 0
             i
         end
    end

end