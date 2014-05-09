% recomputes coverage during part selection
% we don'r re-order coverage withing this function

% X: left center right
% statistics: center, left, .. , right, .. 

function [coverage] = recomputeCoverage(statistics, abstractionTable, outputCoords, coverage, X, startX, recompute, dx, dy, list_depth, lenF, n2Clusters)

    disp('recomputing parts coverage...');
    
    coverageTempEls = []; % format is: [el, cover]
    coverageTempCov = [];
    
    lenStat = size(statistics, 1);
    lenCombs = length(coverage);
    coverage(startX:end) = coverage(startX:end) * 0;
    
    triples = zeros(n2Clusters,n2Clusters,n2Clusters);
   
    is_abstraction = false;
    if ~isempty(abstractionTable)
        is_abstraction = true;
    end
    
    % define parts which need to be re-computed: indsR
    indsR = [];
    indsS = [];
    inds3 = [2,1,4]; % to extract columns of statistics
    statistics = statistics(:, inds3);
    
    if sum(recompute) == n2Clusters
        allRequired = true;     % all parts require re-computation
    else 
        allRequired = false;
    end
    
    if ~allRequired
        
        % Extract X required for re-computation
        parfor i = startX:lenCombs
            curEl = X(i,:);       
            if recompute(curEl(1)) == 1 || recompute(curEl(2)) == 1 || recompute(curEl(3)) == 1
                % re-computation required
                indsR = [indsR; i];
            end
        end

        % Extract statistics required for re-computation
        X_r = X(indsR,:);
        lenCombs = size(X_r, 1);

        % now extract subsets of 'statistics' and 'outputCoords' required for
        % re-computation
        parfor i = 1:lenStat
            curEl = statistics(i,:);       
            if recompute(curEl(1)) == 1 || recompute(curEl(2)) == 1 || recompute(curEl(3)) == 1
                % re-computation required
                indsS = [indsS; i];
            end
        end
        statistics = statistics(indsS,:);
        outputCoords = outputCoords(indsS,:);
        lenStat = size(statistics, 1);
        
    else
        X_r = X;
    end
    
    % fill a table triples
    
    for i = 1:lenCombs
        triples(X_r(i, 1), X_r(i, 2), X_r(i, 3)) = i;
    end
    
    % create a table of format [el, Im, x,y]
    outputCoords = [zeros(lenStat, 1), outputCoords];
    
    for i = 1:lenStat
        
       curEl = statistics(i,:);
       
       % look for the same element in the table X and write the elId
       outputCoords(i,1) = triples(curEl(1), curEl(2), curEl(3));
             
    end
    
    % to exclude already selected parts
    firstCol = outputCoords(:,1);
    indEx = firstCol ~= 0;
    outputCoords = outputCoords(indEx, :);
    
    %----------------------------------------------------------------------
    secondColumn = outputCoords(:,2); % [el, Im, x,y]
    [secondColumn, idx] = sort(secondColumn, 'ascend'); % sort for speed up reasons
    outputCoords = outputCoords(idx, :);
    
    disp('computing coverage through all images...');
    
    parfor i = 1:lenF  % for each image
               
        I = imread(list_depth{i});
        r = size(I,1);
        c = size(I,2);
        
        I = double(I);
        I3 = I(:,:,3);
        I3(I3>0) = 1;  % uncovered area. I3 acts as a mask now
        
        I3 = uint16(I3);
        I2 = zeros(r,c);   % image of zeros
        I2 = uint16(I2);
        
        % extract parts corresponding to this image from outputCoords
        indsI = secondColumn == i;
        els = outputCoords(indsI, :);  % all elements of this image
        lenEls = size(els, 1);
        
        if lenEls > 0   % for some images there are no detected elements
            
            prevEl = els(1,1);

            for j = 1:lenEls

                curEl = els(j,1);

                if curEl ~= prevEl % evaluate coverage
                    I2 = I2.*I3;
                    
                    MM = I2(:);
                    cover = length(MM(MM == 999));
                    
                    coverageTempCov = [coverageTempCov; cover];
                    coverageTempEls = [coverageTempEls; prevEl];
                    
                    
                    I2 = uint16(zeros(r,c)); % apply mask
                    prevEl = curEl;
                end

                % project curEl to the image
                x = els(j, 3);
                y = els(j, 4);
                I2(y-dy:y+dy ,x-dx:x+dx) = ones(2*dy+1, 2*dx+1) * 999;

            end
        end
        
       if mod(i,30) == 0
           i
       end
        
    end

    % now write information from coverageTemp to coverage
    ll = length(coverageTempCov);
    coverageTempCov = double(coverageTempCov);
    
    
    if allRequired
        for i = 1:ll
            coverage(coverageTempEls(i)) = coverage(coverageTempEls(i)) + coverageTempCov(i);
        end
    else
        
        for i = 1:ll
            coverage(indsR(coverageTempEls(i))) = coverage(indsR(coverageTempEls(i))) + coverageTempCov(i);
        end
    end
    
    % now sum coverages up according to the abstraction table
    if is_abstraction
        coverageCopy = coverage;
        for i = 1:lenCombs
            if abstractionTable(i) ~= 0;
                coverage(abstractionTable(i)) = coverage(abstractionTable(i)) + coverageCopy(i);
            end
        end
    end

end
    
    
    
    
    
    



