function [] = AnalyzeVocabLevel(datasetName, vocabLevel, graphLevel, categoryArrIdx, categoryNames, levelItr)
    labelIds = [graphLevel.labelId]';
    imageIds = [graphLevel.imageId]';
    realCount = numel(graphLevel);
    levelIds = repmat(levelItr, realCount, 1);
    exportArr = [labelIds, levelIds, imageIds];
    
    if ~isempty(strfind(datasetName, 'MNIST'))
        categoryNames = cellfun(@(x) str2double(x), categoryNames);
        customOrder = zeros(size(categoryNames,1),1);
        customOrder(categoryNames+1) = 1:size(categoryNames,1);
    else
        customOrder = [];
    end
    
    % Extract simple features from images.
    featureDim = numel(vocabLevel);
    allActivations = unique(exportArr, 'rows');
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
    
    categoryConfMax = zeros(numberOfCategories,numberOfCategories);
    for featureItr = 1:size(catDistArr,1)
        entries = find(catDistArr(featureItr,:));
        if numel(entries)>1
            pairs = combnk(entries, 2);

            for pairItr = 1:size(pairs,1)
                newVal = min(catDistArr(featureItr,pairs(pairItr,:))); 
                categoryConfMax(pairs(pairItr,1), pairs(pairItr,2)) = categoryConfMax(pairs(pairItr,1), pairs(pairItr,2)) + newVal;
                categoryConfMax(pairs(pairItr,2), pairs(pairItr,1)) = categoryConfMax(pairs(pairItr,2), pairs(pairItr,1)) + newVal;
            end
        end

        for entryItr = entries
            categoryConfMax(entryItr, entryItr) = categoryConfMax(entryItr, entryItr) + catDistArr(featureItr, entryItr);
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