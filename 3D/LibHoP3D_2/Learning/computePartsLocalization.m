function [localizationOut] = computePartsLocalization(list_depth, statistics3LayerSieved, statistics3LayerAggregated, nPrevClusters)
    
    
    load(statistics3LayerSieved);  % 'statistics', 'cluster3Depths', 'outputCoords'
    load(statistics3LayerAggregated); % 'X' ,'frequencies', 'curTS', 'triples'
    
    lenF = length(list_depth);
    lenStat = size(statistics, 1);
    lenCombs = size(X, 1);
    
    % define parts which need to be re-computed: indsR
    inds3 = [2,1,4]; % to extract columns of statistics
    statistics = statistics(:, inds3);
    
    triples = zeros(nPrevClusters,nPrevClusters,nPrevClusters);    
    % fill a table triples 
    for i = 1:lenCombs
        triples(X(i, 1), X(i, 2), X(i, 3)) = i;
    end
    
    % create a table of format [el, Im, x,y]
    outputCoords = [zeros(lenStat, 1), outputCoords];
    
    parfor i = 1:lenStat 
       % look for the same element in the table X and write the elId
       outputCoords(i,1) = triples(statistics(i,1), statistics(i,2), statistics(i,3));     
    end
    
    secondColumn = outputCoords(:,2); % [el, Im, x,y]
    [secondColumn, idx] = sort(secondColumn, 'ascend'); % sort for speed up
    outputCoords = outputCoords(idx, :);
    
    disp('computing localization of parts through all images...');
 
    localizationIntermediate = zeros(lenCombs, lenF);
    
    for i = 1:lenF
        I = imread(list_depth{i});
        r = size(I,1);
        c = size(I,2);
        
        I3 = zeros(r,c);
        
        localizationLocal = zeros(lenCombs, 1);
        numElements = zeros(lenCombs, 1);
        
        % extract parts corresponding to this image from outputCoords
        indsI = secondColumn == i;
        els = outputCoords(indsI, :);
        
        firstColumn = els(:, 1);
        [~, inddds] = sort(firstColumn);
        els = els(inddds, :);
        lenEls = size(els, 1);
        
        if lenEls > 0   % for some images there are no detected elements
            
            % project all elements to the data
            for j = 1:lenEls
                I3(els(j, 4), els(j, 3)) = els(j, 1);
            end
            
            % consider a local neigbourhood of each element
            for j = 1:lenEls
                curEl = els(j, 1);
                y = els(j, 4);
                x = els(j, 3);
                window = I3(y-2 : y+2, x-2 : x+2);
                
                a = window(window == curEl);
                count = length(a);
                localizationLocal(curEl) = localizationLocal(curEl) + count;
                numElements(curEl) = numElements(curEl) + 1;
            end
            
        end
        
            for j = 1:lenCombs
                if numElements(j) > 0
                    localizationLocal(j) = localizationLocal(j)/ numElements(j);
                end
            end
        localizationIntermediate(:, i) = localizationLocal;
    end
    
    localizationOut = zeros(lenCombs, 1);
    for j = 1:lenCombs
        curLine = localizationIntermediate(j, :);
        localizationOut(j) = sum(curLine) / nnz(curLine);  % average
    end

end

