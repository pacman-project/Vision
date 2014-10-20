function [] = AnalyzeNodes(datasetName)

    % Load relevant info.
    load([pwd '/output/' datasetName '/vb.mat']);
    load([pwd '/output/' datasetName '/export.mat']);
    if ~isempty(strfind(datasetName, 'MNIST'))
        categoryNames = cellfun(@(x) str2double(x), categoryNames);
        customOrder = zeros(size(categoryNames,1),1);
        customOrder(categoryNames+1) = 1:size(categoryNames,1);
    else
        customOrder = [];
    end

    % Extract simple features from images.
    featureDims = cellfun(@(x) numel(x), vocabulary);
    featureDim = sum(featureDims);
    cumSums = int32(cumsum(featureDims));
    cumSums(2:end) = cumSums(1:(end-1));
    cumSums(1) = 0;
    cumSums(end+1) = featureDim;
    allActivations = exportArr(:,[1,4,5]);
    allActivations(:,1) = allActivations(:,1) + cumSums(allActivations(:,2));
    allActivations = unique(allActivations, 'rows');
    allActivations(:,4) = int32(categoryArrIdx(allActivations(:,3)));
    
    numberOfCategories = max(categoryArrIdx);
    categoryArrList = 1:max(categoryArrIdx);
    catDistArr = zeros(featureDim, size(categoryArrList,2));
    % Analyze each feature's discriminative properties.
    for featureItr = 1:featureDim
        activations = allActivations(allActivations(:,1) == featureItr,:);
        catDist = hist(activations(:,4), categoryArrList);
        catDistArr(featureItr,:) = catDist;
    end
    
    for levelItr = 1:max(allActivations(:,2))
        categoryConfMax = zeros(numberOfCategories,numberOfCategories);
        levelDistArr = catDistArr((cumSums(levelItr)+1):cumSums(levelItr+1),:); 
        for featureItr = 1:size(levelDistArr,1)
            entries = find(levelDistArr(featureItr,:));
            if numel(entries)>1
                pairs = combnk(entries, 2);

                for pairItr = 1:size(pairs,1)
                    newVal = min(levelDistArr(featureItr,pairs(pairItr,:))); 
                    categoryConfMax(pairs(pairItr,1), pairs(pairItr,2)) = categoryConfMax(pairs(pairItr,1), pairs(pairItr,2)) + newVal;
                    categoryConfMax(pairs(pairItr,2), pairs(pairItr,1)) = categoryConfMax(pairs(pairItr,2), pairs(pairItr,1)) + newVal;
                end
            end

            for entryItr = entries
                categoryConfMax(entryItr, entryItr) = categoryConfMax(entryItr, entryItr) + levelDistArr(featureItr, entryItr);
            end
        end
        if ~isempty(customOrder)
            categoryConfMax = categoryConfMax(customOrder, customOrder);
        end
        if ~exist([pwd '/categorization/analysis/' datasetName], 'dir')
            mkdir([pwd '/categorization/analysis/' datasetName]);
        end
        save([pwd '/categorization/analysis/' datasetName '/level' num2str(levelItr) '.mat'], 'categoryConfMax');
        imwrite(categoryConfMax/max(max(categoryConfMax)), [pwd '/categorization/analysis/' datasetName '/level' num2str(levelItr) '.png']);
    end
end