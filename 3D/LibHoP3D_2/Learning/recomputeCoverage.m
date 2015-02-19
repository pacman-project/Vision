% recomputes coverage during part selection
% we don'r re-order coverage withing this function

% X: left center right
% statistics: center, left, .. , right, .. 

% numRecompute - shows how many elements we need to recompute

function [coverage] = recomputeCoverage(statistics, outputCoords, coverage, X, startX, dx, dy,...
                                        list_depth, lenF, nPrevClusters, numRecompute, is_sparse, triples)

    disp('recomputing parts coverage...');
    
    X = uint32(X);
    statistics = uint32(statistics);
    
    lenStat = size(statistics, 1);
    lenCombs = length(coverage);
    
    if startX+numRecompute-1 > lenCombs
        numRecompute = lenCombs - startX + 1;
    end
        
        
    coverage(startX:startX+numRecompute-1) = coverage(startX:startX+numRecompute-1) * 0;
    
    
    % create a matrix
    if isempty(triples)
        triples = zeros(nPrevClusters, nPrevClusters, nPrevClusters);
    end
    
    % define parts which need to be re-computed: indsR
    inds3 = [2,1,4]; % to extract columns of statistics
    statistics = statistics(:, inds3);
        
    
    
    % create a table of format [el, Im, x,y]
    outputCoords = [zeros(lenStat, 1), outputCoords];
    
    if is_sparse 
        % fill a table triples
        for i = startX:startX+numRecompute-1  % lenCombs
            triples{X(i, 1)}(X(i,2), X(i,3)) = i;
        end
        
        % fill a table outputCoords
        for i = 1:lenStat
            outputCoords(i,1) = triples{statistics(i,1)}(statistics(i,2), statistics(i,3)); 
        end
               
    else    % if it is NOT sparse
        
        % fill a table triples
        inds = startX:1:startX+numRecompute-1;
        indSS = sub2ind(size(triples), X(inds,1), X(inds,2), X(inds,3)); 
        triples(indSS) = inds;
    
%         for i = startX:startX+numRecompute-1  % lenCombs
%             triples(X(i, 1), X(i,2), X(i,3)) = i;
%         end
        
        % fill a table outputCoords
        
        inds = 1:1:lenStat;
        indSS = sub2ind(size(triples), statistics(inds,1), statistics(inds,2), statistics(inds,3)); 
        outputCoords(inds, 1) = triples(indSS);
        
%         for i = 1:lenStat 
%             outputCoords(i,1) = triples(statistics(i,1), statistics(i,2), statistics(i,3));  
%         end 
    end

    % to exclude already selected parts (and those which do not need to be recomputed)
    firstCol = outputCoords(:,1);
    indEx = firstCol ~= 0;
    outputCoords = outputCoords(indEx, :);
    
    %----------------------------------------------------------------------
    secondColumn = outputCoords(:,2); % [el, Im, x,y]
    [secondColumn, idx] = sort(secondColumn, 'ascend'); % sort for speed up
    outputCoords = outputCoords(idx, :);
    
    disp('computing coverage through all images...');
    sumCov = zeros(numRecompute, lenF);

    parfor i = 1:lenF  % for each image
               
        I = imread(list_depth{i});
        r = size(I,1);
        c = size(I,2);
        
        I3 = I(:,:,3);
        I3(I3>0) = 1;  % uncovered area. I3 acts as an object mask now
        
        I3 = uint8(I3);
        I2 = uint8(zeros(r,c));   % image of zeros
        localSumCov = zeros(numRecompute, 1);
        
        % extract parts corresponding to this image from outputCoords
        indsI = secondColumn == i;
        els = outputCoords(indsI, :);  % all elements of this image
        firstColumn = els(:, 1);
        [~, inddds] = sort(firstColumn);
        els = els(inddds, :);
        lenEls = size(els, 1);
        
        if lenEls > 0   % for some images there are no detected elements   
            prevEl = els(1,1);
            
            for j = 1:lenEls
                
                curEl = els(j,1);

                if curEl ~= prevEl % evaluate coverage of the prevEl
                    
                    I2 = I2+I3;
                    A = I2 == 2;
                    cover  = nnz(A);
                    
                    localSumCov(prevEl - startX + 1) = localSumCov(prevEl - startX + 1) + cover;
                    
                    I2 = uint8(zeros(r,c)); % apply mask
                    prevEl = curEl;
                end

                % project curEl to the image
                x = els(j, 3);
                y = els(j, 4);
                I2(y-dy:y+dy ,x-dx:x+dx) = ones(2*dy+1, 2*dx+1);
    
            end
        end
        
        sumCov(:, i) = localSumCov;
        
        if mod(i,30) == 0
           i
        end  
    end
    
    % now write information from coverageTemp to coverage
    coverage(startX:startX+numRecompute-1) = sum(sumCov,2);
end






%     lineLength = 5;
%     
%     
%     parfor i = 1:lenF/lineLength  % for each image
%         
%         imInds = 1+(i-1)*lineLength:1:i*lineLength;
%         sizeC = zeros(1,lineLength);
%         sizeR = zeros(1,lineLength);
%         I = {};
%         
%         for j = 1:lineLength
%             I{j} = imread(list_depth{imInds(j)});
%             sizeR(j) = size(I{j}, 1);
%             sizeC(j) = size(I{j}, 2);           
%         end
%         
%         % concatenate matrices
%         rBig = sum(sizeR);
%         cBig = max(sizeC);
%         
%         offsets = [0, cumsum(sizeR(1:4))];
% 
%         I3 = gpuArray.false(rBig, cBig);
%         I2 = gpuArray.false(rBig, cBig);   % image of zeros
%         for j = 1:lineLength
%             curI = I{j};
%             curI3 = curI(:,:,3);
%             curI3(curI3>0) = 1;
%             curI3 = logical(curI3);
%             I3(offsets(j) + 1: offsets(j) + sizeR(j), 1: sizeC(j)) = curI3;
%         end
% %         I33 = gather(I3);
% %         imshow(I33);
%    
%         % extract parts corresponding to this image from outputCoords
%         indsI = (secondColumn >= imInds(1)) & (secondColumn <= imInds(lineLength));
%         els = outputCoords(indsI, :);  % all elements of these images
%         for ii = 2:lineLength
%             inds = els(:,2) == ii;
%             els(inds, 4) = els(inds, 4) + offsets(ii);
%         end
%         
%         firstColumn = els(:, 1);
%         [~, inddds] = sort(firstColumn);
%         els = els(inddds, :);
%         lenEls = size(els, 1);
%         
%         if lenEls > 0   % for some images there are no detected elements
%             
%             prevEl = els(1,1);
% 
%             for j = 1:lenEls
%                 
%                 curEl = els(j,1);
% 
%                 if curEl ~= prevEl % evaluate coverage of the prevEl
%                     
%                     
% %                     I2 = I2 & I3;
% %                     
% %                     
% %                     A = I2 == 9999;
% %                     cover  = nnz(A);
% % %                     
% %                     coverageTempCov = [coverageTempCov; cover];
% %                     coverageTempEls = [coverageTempEls; prevEl];
%                     
%                     I2 = gpuArray.false(rBig, cBig); % apply mask
%                     prevEl = curEl;
%                 end
% 
%                 % project curEl to the image
%                 x = els(j, 3);
%                 y = els(j, 4);
%                 I2(y-dy:y+dy ,x-dx:x+dx) = true(2*dy+1, 2*dx+1);
%                 
%             end
%         end
%         
%         if mod(i,10) == 0
%            i
%         end
%         
%     end
%     
% 
%     % now write information from coverageTemp to coverage
%     ll = length(coverageTempCov);
%     coverageTempCov = double(coverageTempCov);
%     
%     for i = 1:ll
%         coverage(coverageTempEls(i)) = coverage(coverageTempEls(i)) + coverageTempCov(i);
%     end
% end



    
    
    
    
    
    



