% this is to compute coverage of all vertical layers (3,5, etc)

% X - are combinations 
function [] = computeCoverage3_4(list_El, lenF, nClusters, n2Clusters, n3Clusters, X, partsOut, displacement, abstractionLevel,  ...
               areLayersRequired, outRoot, combs, largestLine, wCoverage, wOverlap, isInhibitionRequired, is_downsampling, dowsample_rate, elPath)

   
    is_third_layer = areLayersRequired(3);
    is_4th_layer = areLayersRequired(4);
    is_inhibition_3_layer = isInhibitionRequired(3);
    is_inhibition_4_layer = isInhibitionRequired(4); 
    
    hDispl = round(displacement/2); % half displ
    
    lenDPW = length(elPath);
    
%     isTrim = true;
%     isErrosion = false;
%     isX = false;
%     isY = false;
    
    
    table = zeros(n2Clusters, 2);
    for i = 1:n2Clusters
        [clusterX, clusterY] = compute2derivatives(i, nClusters);
        table(i, 1) = clusterX;
        table(i, 2) = clusterY;
    end
    
    table3 = uint16(zeros(n2Clusters,n2Clusters,n2Clusters));
    
    % compute the the abstraction table
    abstractionTable = computeAbstractionTable(X, nClusters, n2Clusters, abstractionLevel);
    
    lenAbst = length(abstractionTable);
    
    for i = 1:n3Clusters
        cur = partsOut(i,:);  % lest centre right
        table3(cur(1),cur(2),cur(3)) = i;   % all selected parts are now in this table
    end
    
    table3Abst = table3;  % assume elements refer to themselve
    
    for i = 1:lenAbst      
        if abstractionTable(i) > 0  % part X(i,:) refers to part X(pointer, :)
            pointer = abstractionTable(i);
            curX = X(pointer,:); % element we refer to

            table3Abst(X(1), X(2), X(3)) = table3(curX(1),curX(2),curX(3));           
        end
    end
    
    
   parfor i = 1:lenF 
        
%         I = imread(list_depth{i}); 
%         I = double(I);
%         if is_downsampling
%             I = imresize(I, dowsample_rate);
%         end
        marks2 = imread(list_El{i});
        
        % trim image I
%         [I, ~, ~, ~, r, c, is_successfull] = preliminaryProcessing(I, isErrosion, discSize, isX, isY, isTrim, dxKernel, sigmaKernelSize, sigma);
        
        [r, c] = size(marks2);
%         if r~=rEl || c ~= cEl
%             disp('ERROR');
%         end
        
        marks3 = zeros(r,c);
        % Imarks = zeros(r,c,3);
        
        [rows, cols] = find(marks2 > 0);
        nEl = length(rows);
        
        for j = 1: nEl
            
            % check whether it is close to the boundary
            if rows(j) < displacement + 1  || rows(j) > r - displacement  || cols(j) < displacement + 1 || cols(j) > c - displacement
                continue;
            else
                % otherwise try to match something around this object
                
                central = marks2(rows(j), cols(j));
%                 depth_centr =  I(rows(j), cols(j));
                
                
                % check what are left and right neighbours
                
                lefts =  [marks2(rows(j), cols(j) - displacement), marks2(rows(j), cols(j) - hDispl)];  % , marks2(rows(j), cols(j) - displ - hDispl)
                rights = [marks2(rows(j), cols(j) + displacement), marks2(rows(j), cols(j) + hDispl)];  % , marks2(rows(j), cols(j) + displ + hDispl)
                indsL = find(lefts > 0);
                indsR = find(rights > 0);

                if (isempty(indsL)) || (isempty(indsR))  % nothing can be matched
                    continue;
                end
%                 depthsLeft =  [I(rows(j), cols(j) - displ), I(rows(j), cols(j) - hDispl)];
%                 depthsRight = [I(rows(j), cols(j) + displ), I(rows(j), cols(j) + hDispl)];
                lefts = lefts(indsL);
                rights = rights(indsR);
                
%                 depthsLeft = depthsLeft(indsL);
%                 depthsRight = depthsRight(indsR);
                
                done = false; % local structure if matched to the vocabulary elements
                ii = 1;
                jj = 1;
                
                while (~done && ii <= length(indsL) && jj <= length(indsR))
                    left = lefts(ii);
                    right = rights(jj);
%                     dLeft = depthsLeft(ii);
%                     dRight = depthsRight(jj);

                    
                    el = [left, central, right];
                    
                    % matching
                    
                    curEl = table3Abst(el(1), el(2), el(3));  % all or nodes are already in this table
                    
                    if curEl ~= 0
                        marks3(rows(j), cols(j)) = curEl;
                        done = true;                     
                    end
                    
                    jj = jj+1; % increment loop variable
                    if jj > length(indsR)
                        jj = 1;
                        ii = ii+1;
                    end
                    
                end
                               
            end 
        end
        % save the image
        
        curStr = list_El{i};

        fileName = curStr(lenDPW+1:end);
        outFile = [outRoot, fileName];

        ll = strfind(outFile, '/');
        ll = ll(end); % last position
        folderName = outFile(1:ll);
        b = exist(folderName,'dir');

        if b == 0
            mkdir(folderName);
        end

        marks3 = uint16(marks3);
        imwrite(marks3, outFile, 'png');
        
        if mod(i,200) == 0
            i
        end
        
   end
        
        
%         % ----------------------RECOGNITION OF THE 4TH LAYER HERE----------

  % table4 = uint16(zeros(n3Clusters, n3Clusters, n3Clusters));
%           marks4 = zeros(r,c);
%         
%         [rows, cols] = find(marks3 > 0);
%         nEl = length(rows);

%     for iii = 1:n4Clusters
%         cur = triples4Out(iii,:);
%         table4(cur(1),cur(2),cur(3)) = iii;
%     end
%         
%         
%         for j = 1:nEl
%             
%             % check whether it is close to the boundary
%             if rows(j) < displ + 1  || rows(j) > r - displ  || cols(j) < displ + 1 || cols(j) > c - displ
%                 continue;
%             else
%                 central = marks3(rows(j), cols(j));
% %                 depth_centr =  I(rows(j), cols(j));
%                 
%                 
%                 % check what are top and bottom neighbours
%                 
%                 tops =  [marks3(rows(j) - displ, cols(j)), marks3(rows(j) - hDispl, cols(j)), ...
%                          marks3(rows(j) - displ, cols(j) - hDispl), marks3(rows(j) - displ, cols(j) + hDispl) ]; 
%                 
%                 bottoms = [marks3(rows(j) + displ, cols(j)),  marks3(rows(j) + hDispl, cols(j)), ...
%                          marks3(rows(j) + displ, cols(j) - hDispl), marks3(rows(j) + displ, cols(j) + hDispl)];
%                 
% %                 tops =  [marks3(rows(j) - displ, cols(j)),          marks3(rows(j) - hDispl, cols(j)), ...
% %                          marks3(rows(j) - displ, cols(j) - hDispl), marks3(rows(j) - displ, cols(j) + hDispl), ...
% %                          marks3(rows(j) - hDispl, cols(j)- hDispl), marks3(rows(j) - hDispl, cols(j)+ hDispl)]; 
% %                 
% %                 bottoms = [marks3(rows(j) + displ, cols(j)),          marks3(rows(j) + hDispl, cols(j)), ...
% %                          marks3(rows(j) + displ, cols(j) - hDispl), marks3(rows(j) + displ, cols(j) + hDispl), ...
% %                          marks3(rows(j) + hDispl, cols(j)- hDispl), marks3(rows(j) + hDispl, cols(j)+ hDispl)];
%                      
%                      
%                 indsT = find(tops > 0);
%                 indsB = find(bottoms > 0);
% 
%                 if (isempty(indsT)) || (isempty(indsB))
%                     continue;
%                 end
%                 
% %                 depthsTop =    [I(rows(j) - displ, cols(j)), I(rows(j) - hDispl, cols(j))];
% %                 depthsBottom = [I(rows(j) + displ, cols(j)), I(rows(j) + hDispl, cols(j))];
% 
%                 tops = tops(indsT);
%                 bottoms = bottoms(indsB);
%                 
% %                 depthsTop = depthsTop(indsT);
% %                 depthsBottom = depthsBottom(indsB);
%                 
%                 done = false; % local structure if matched to the vocabulary elements
%                 ii = 1;
%                 jj = 1;
%                 
%                 while (~done && ii <= length(indsT) && jj <= length(indsB))
%                     top = tops(ii);
%                     bottom = bottoms(jj);
% %                     dTop = depthsTop(ii);
% %                     dBottom = depthsBottom(jj);
%                     
%                     el = [top, central, bottom];
%                     
% %                     % this is imitation of the OR-nodes
% %                     el = orNode4(el, table3, triples3Out);
%                     
%                     % matching
%                     curEl = table4(el(1), el(2), el(3));
%                     
%                     if curEl ~= 0
%                         marks4(rows(j), cols(j)) = curEl;
%                         done = true;                 
%                     end
%                     
%                     jj = jj+1; % increment loop variable
%                     if jj > length(indsB)
%                         jj = 1;
%                         ii = ii+1;
%                     end
%                     
%                 end
%                                
%             end 
%         end
%         
%         %---------------------------end of 4th layer recognition-----------
%         
%         if is_third_layer && is_4th_layer  % write the results somewhere
%             
%             shift = n2Clusters;
%             marks3 = marks3 + shift;
%             marks3(marks3 == shift) = 0;
%             
% %             imtool(I, [min(min(I)), max(max(I))]);
% %             imtool(marks3, [min(min(marks3)), max(max(marks3))]);
% %             imtool(marks4, [min(min(marks4)), max(max(marks4))]);
%             
% %             len = length(marks3>0)
%             
%             shift = n2Clusters + n3Clusters;
%             marks4 = marks4 + shift;
%             marks4(marks4 == shift) = 0;
%             
%             curStr = list_El{i};
%     %       ll = strfind(curStr, '/');
%     %       fileName = curStr(ll:end);
% 
%             fileName = curStr(lenDPW+1:end);
%             outFile = [outRoot, fileName];
% 
%             ll = strfind(outFile, '/');
%             ll = ll(end); % last position
%             folderName = outFile(1:ll);
%             b = exist(folderName,'dir');
% 
%             if b == 0
%                 mkdir(folderName);
%             end
%             
%             % create many dimensional feature vector
%             Imarks(:,:,1) = marks2;
%             Imarks(:,:,2) = marks3;
%             Imarks(:,:,3) = marks4;
%             Imarks = uint16(Imarks);
%             
%             imwrite(Imarks, outFile, 'png');
%         end
% 
%       %  imtool(marks3, [min(min(marks3)), max(max(marks3))]);
%         if mod(i,10) == 0
%             i
%         end
%     end
%     

        
end