% this function computes measure of part's localization in surfaces

% the formula for localization is the simplest one could imagine:
% localizationOut = (number of occurrences of part realizations)/(number of occurrences in neighbourhoods of these part realizations);

function [localizationOut] = computePartsLocalization(list_depth, statistics3LayerSieved, statistics3LayerAggregated, nPrevClusters, winSize)
    
    lenF = length(list_depth);
    
    load(statistics3LayerSieved);  % 'statistics', 'cluster3Depths', 'outputCoords'

    
    b = load(statistics3LayerAggregated); % 'X' ,'frequencies', 'curTS', 'triples'
    X = b.X;
    
    clear('a', 'b');
    
    lenStat = size(statistics, 1);
    lenCombs = size(X, 1);
    
    if nPrevClusters > 700
        is_sparse = true;
    else
        is_sparse = false;
    end
    
    % define parts which need to be re-computed: indsR
    inds3 = [2,1,4]; % to extract columns of statistics
    
    statistics = statistics(:, inds3);   % do not consider depths
    
    
    if is_sparse
        parfor i = 1:nPrevClusters
            triples{i} = sparse(zeros(nPrevClusters, nPrevClusters));
        end
    else
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
    
    % create a table of format [el, Im, x,y]
    outputCoords = [zeros(lenStat, 1), outputCoords];
    
    if ~is_sparse
        parfor i = 1:lenStat 
           % look for the same element in the table X and write the elId
           outputCoords(i,1) = triples(statistics(i,1), statistics(i,2), statistics(i,3)); 
           
           if mod(i, 100000) == 0
               i
           end
        end
    else
        parfor i = 1:lenStat 
           % look for the same element in the table X and write the elId
           outputCoords(i,1) = triples{statistics(i,1)}(statistics(i,2), statistics(i,3)); 
           if mod(i, 100000) == 0
               i
           end
        end  
    end
    
    secondColumn = outputCoords(:,2); % [el, Im, x,y]
    [secondColumn, idx] = sort(secondColumn, 'ascend'); % sort according to image id
    outputCoords = outputCoords(idx, :);
    
    disp('computing localization of parts through all images...');
    
    localizationGlobal = zeros(lenCombs, 1);
    numElementsGlobal = zeros(lenCombs, 1);
    
    for i = 1:lenF
        
        I = imread(list_depth{i});
        r = size(I,1);
        c = size(I,2);
        
        I3 = zeros(r,c);
        
        localizationLocal = zeros(lenCombs, 1);
        numElementsLocal= zeros(lenCombs, 1);
    
        % extract parts corresponding to this image from outputCoords
        indsI = secondColumn == i;
        els = outputCoords(indsI, :);    % should be: el, image, x,y.
        
        firstColumn = els(:, 1);
        [~, inddds] = sort(firstColumn);
        els = els(inddds, :);
        lenEls = size(els, 1);  % all the elements corresponding to this image
        
        if lenEls > 0   % for some images there are no detected elements
            
            % project all elements to the data
            for j = 1:lenEls
                I3(els(j, 4), els(j, 3)) = els(j, 1);
            end
            
            % consider a local neigbourhood of each part
            for j = 1:lenEls
                curEl = els(j, 1);
                y = els(j, 4);
                x = els(j, 3);
                % check image borders
                topBorder = y-winSize;
                bottomBorder = y+winSize;
                leftBorder = x-winSize;
                rightBorder = x+winSize;
                
                if topBorder <= 0
                    topBorder = 1;
                end
                if leftBorder <= 0
                    leftBorder = 1;
                end
                if bottomBorder > r
                    bottomBorder = r;
                end
                if rightBorder > c
                    rightBorder = c;
                end
                
                    
                window = I3(topBorder : bottomBorder, leftBorder : rightBorder);
                
                a = window(window == curEl);
                count = length(a);
                localizationLocal(curEl) = localizationLocal(curEl) + count;
                numElementsLocal(curEl) = numElementsLocal(curEl) + 1;
            end
          
            
        end
        
        if mod(i, 10) == 0
            i
        end
        
        localizationGlobal = localizationGlobal + localizationLocal;
        numElementsGlobal = numElementsGlobal + numElementsLocal;
    end
    
    localizationOut = numElementsGlobal ./ (localizationGlobal + numElementsGlobal);  

end

