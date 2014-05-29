function finalImg = GetActivationImage(exportArr, vocabCounts)
    global vocabulary;
    %% Initialize variables.
    imageSize = [480, 640];
    partIds = exportArr(:,1);
    levelIds = exportArr(:,4);
    numberOfDetectedLevels = max(levelIds);
    compSize = 5;
    sidePadding = 40;

    %% Get the ratio of compositions to be preserved.
    maxCompsPerLevel = floor(imageSize(:,2)/(compSize+1));
    preserveRatio = maxCompsPerLevel/max(vocabCounts);
    numberOfLevels = numel(vocabCounts);
    
    %% Create level-wise images.
    activationImg = zeros(imageSize);
    previousVocabLevelPos = [];
    prevRealIdArr = [];
    for levelItr = 1:numberOfDetectedLevels
        preservedComps = vocabCounts(levelItr);
        if levelItr ~= 1
           preservedComps = floor(preservedComps * preserveRatio);
        end
        
        if preservedComps == 0
            break;
        end
        
        activatedParts = unique(partIds(levelIds==levelItr & partIds <= preservedComps));
        
        levelImg = zeros(compSize, preservedComps * (compSize + 1));
        isActivated = ismember(1:preservedComps, activatedParts);
        currentVocabLevelPos = zeros(preservedComps, 2);
        for compItr = 1:preservedComps
            if isActivated(compItr)
                value = 2;
            else
                value = 1;
            end
            beginOffset = ((compItr-1) * (compSize + 1)) + 1;
            levelImg(:,beginOffset:(beginOffset + compSize - 1)) = value;
            currentVocabLevelPos(compItr,:) = [round(compSize/2)-1, beginOffset+round(compSize/2)-1];
            
            realIdArr = 1:numel(isActivated);
            realIdArr(isActivated) = 1:nnz(isActivated);
            realIdArr(~isActivated) = 0;
        end
        
        % If no nodes have been activated, exit.
        if nnz(isActivated) == 0
            break;
        end
        
        % Place the level image to its right spot.
        emptySpace = floor(((imageSize(:,1) - sidePadding*2) - (numberOfLevels*compSize))/(numberOfLevels-1));
        offset = (levelItr-1)*(emptySpace + compSize) + 1 + sidePadding;
        leftOffset = floor((imageSize(:,2)-size(levelImg,2))/2)+1;
        activationImg((end-offset):(((end-offset) + size(levelImg,1))- 1), leftOffset:((leftOffset+size(levelImg,2))-1)) = levelImg;
        currentVocabLevelPos = currentVocabLevelPos + repmat([imageSize(1)-offset, leftOffset], preservedComps, 1);
        
        % eliminate unused current comps
        currentVocabLevelPos = currentVocabLevelPos(isActivated',:);
        
        %% Create edges between the layers
        currentVocabLevelChildren = {vocabulary{levelItr}(1:preservedComps).children};
        if levelItr>1
            prevPartIds = prevRealIdArr(prevRealIdArr>0);
            
            % If no nodes are activated in the previous layer, stop.
            if isempty(prevPartIds)
                break;
            end
            
            currentVocabLevelChildren = currentVocabLevelChildren(isActivated);
            
            % Eliminate compositions with large ids.
            currentVocabLevelChildren = cellfun(@(x) x(ismember(x, prevPartIds)), currentVocabLevelChildren, 'UniformOutput', false);
            edges = cell(numel(currentVocabLevelChildren),1);
            for compItr = 1:numel(currentVocabLevelChildren)
                children = currentVocabLevelChildren{compItr};
                
                % If no valid children exist, move on.
                if isempty(children)
                    continue;
                end
                
                % Link the nodes.
                edges{compItr} = [repmat(currentVocabLevelPos(compItr,:), numel(children),1), ...
                    previousVocabLevelPos(children,:)];
            end
            edges = cat(1, edges{:});
            
            % Create lines between the nodes.
            if isempty(edges)
                break;
            end
            ind = drawline(edges(:,1:2), edges(:,3:4), imageSize);
            [xInd, ~] = ind2sub(imageSize, ind');
            
            % Eliminate parts of lines too close to the both levels.
            validInd = xInd<(previousVocabLevelPos(1,1)-2) & xInd>(currentVocabLevelPos(1,1)+2);
            ind = ind(validInd);
            activationImg(ind) = 3;
        end
        previousVocabLevelPos = currentVocabLevelPos;
        prevRealIdArr = realIdArr;
    end
    
    % Give the activation image some colored outlook.
    finalImg = label2rgb(activationImg, [0.8,0,0; 0, 1, 0; 0, 1, 0], 'k');
end

