% this is the main function to co-occurrence statistics of corners and
% edgese


% output is an array of type:
% [elCentral,  el,z,   el,z,  el,z,   el,z,   el,z,   el,z,  el,z,   el,z ] - elements and relative depths
% z - are relative positions w.r.t the central element
% isFull is true when all elements exist in a row
% fieldSize is [sizeX, sizeY, sizeZ]

function [] = CollectStats_NextLayersCorners(list_els, list_input, list_mask, lenF, filtOptions,  ...
                                                        nPrevClusters, layerID, depthStep, dataSetNumber, inputDataType, zScale, ...
                                                        maxRelDepth, cluster1Bounds, nClusters, options, partIDStart, inPath, outRoot)
    
    disp('collecting statistics of corners and edges...');                                                
                                                    
    is_visualization = true;
    structureTensorRad = 12;
    is2D = false;
    nonVanishThresh = 0.1;
    
    if dataSetNumber ~= 2
        list_mask = zeros(1, lenF);
    end

    lenDPW = length(inPath);
    
    parfor i = 1:lenF % This is done to speed up                    

            I = imread(list_input{i});     % depth image
            marksPrev = imread(list_els{i});  % elements of the previous layer
            % read image with normals
            normFile = [list_els{i}(1:end - 4), '_N.mat']; 
            Normals = load(normFile); 
            Normals = Normals.Normals;
            
            Nx = Normals(:,:,1);
            Ny = Normals(:,:,2);
            Nz = Normals(:,:,3);
            
            %% outputFile
            
            curStr = list_input{i};
            fileName = curStr(lenDPW+1:end);
            outFile = [outRoot, fileName];

            ll = strfind(outFile, '/');
            lll = strfind(outFile, '\');
            ll = [ll, lll];
            ll = max(ll); % last position
            folderName = outFile(1:ll);
            b = exist(folderName,'dir');

            if b == 0
                mkdir(folderName);
            end

            % chech if file exists
            b = exist(outFile, 'file');

%              if b
%                  continue;
%              end
            
            I = I(:,:,1);
            I = I * zScale;
            [r,c] = size(I);
            Iout = zeros(r,c,3);
            markOut = zeros(r,c);
            
            if dataSetNumber == 2
                mask = imread(list_mask{i});
            else
                mask = [];
            end
            
            if dataSetNumber == 1 || dataSetNumber == 3 || dataSetNumber == 4 % Aim@Shape dataset || Vladislav_STD 
                [I, ~, ~, mask, is_successfull] = preliminaryProcessing(I, mask, filtOptions);                
            elseif dataSetNumber == 2  % Washington data set 
                [I, ~, ~, mask, is_successfull] = preliminaryProcessing(I, mask, filtOptions);
            end
            if ~is_successfull
                continue;
            end

            [rEl, cEl] = size(mask);  % all three images should be of the same size!
            if r~=rEl || c ~= cEl
                disp('ERROR');
            end
            [rEl, cEl] = size(marksPrev);
            if r~=rEl || c ~= cEl
                disp('ERROR');
            end
            if ~is_successfull
                disp('ERROR');
            end
             

            for x = structureTensorRad+1 : options.step: c - structureTensorRad - 1
                for y = structureTensorRad+1 : options.step: r - structureTensorRad - 1
                    
                    if mask(y,x) == 0
                        continue;
                    end

                    [marksSel, indsX, indsY, ~, ~] = GetPartsNeighbour(I, marksPrev, x, y, I(y,x), structureTensorRad, is2D);
                     
                    if length(indsX) < 3
                        continue;
                    end
                    
                    % get normals at these points
                    inds = sub2ind([r,c], indsY, indsX);
                    
%                     imtool(Normals);
                    
                    Ns_x = Nx(inds)';
                    Ns_y = Ny(inds)';
                    Ns_z = Nz(inds)';
                    Ns = [Ns_x, Ns_y, Ns_z];
                    
                    if std(Ns_x) < 0.1 && std(Ns_y) < 0.1
                        continue;
                    end

                    [T, V, D] = structureTensor3(Ns);
                    b = diag(D);

                    for j =1:3
                        Iout(y,x,j) = b(j); 
                    end
                    
                    if (b(2) >= nonVanishThresh) && b(1) <= nonVanishThresh  % this is an edge 
                        
                        % split all vectors into two clusters
                        
                        try
                            IDs = kmeans(Ns, 2);   
                        catch
                            continue
                        end
                        id1 = IDs == 1;
                        id2 = IDs == 2;
                        Vector1 = mean(Ns(id1, :),1);
                        Vector2 = mean(Ns(id2, :),1);
                        
                        angle = angleVec(Vector1, Vector2);
                        angle = 90 - angle; 
                        curPartID = partIDStart + define1Cluster(angle, cluster1Bounds, nClusters);
                        
                        markOut(y,x) = curPartID;
                        
                    elseif (b(2) >= nonVanishThresh) && b(1) >= nonVanishThresh
                        
                        % split all vectors into tree clusters
                        try
                            IDs = kmeans(Ns, 3);
                        catch
                            a = 2;
                        end
                        id1 = IDs == 1;
                        id2 = IDs == 2;
                        id3 = IDs == 3;
                        Vector1 = mean(Ns(id1, :),1);
                        Vector2 = mean(Ns(id2, :),1);
                        Vector3 = mean(Ns(id3, :),1);
                        
                        angle12 = 90 - angleVec(Vector1, Vector2);
                        angle13 = 90 - angleVec(Vector1, Vector3);
                        angle23 = 90 - angleVec(Vector2, Vector3);
                        
                        curPartID = partIDStart + define1Cluster(angle, cluster1Bounds, nClusters);
                        
                        Iout(y,x) = curPartID;
                        continue;
                        
                    end

                end
            end
            markOut = uint8(markOut);
            imwrite(markOut, outFile, 'png');
%             imtool(markOut, [0, 60]);
%             imtool(Normals);
            
            a = 2;
            
            
   end
                                                   
end

function angleGrad = angleVec(V1, V2) % returns the angle in degrees

    angleGrad = 180 * acos(dot(V1, V2) / (norm(V1) * norm(V2))) / pi;
end

function [af, bf, cf, df] = myFilter4(a,b,c,d, ids)
    af = a(ids);
    bf = b(ids);
    cf = c(ids);
    df = d(ids);
end

function plotFrame(V, centre, vecLen, vectColors, n)
    % plot the eigenvectors in 3D
    for ii = 1:n
        curVect = V(:, ii);
        curColor = vectColors(ii, :);
        XX = [centre(1), centre(1) + vecLen * curVect(1)];
        YY = [centre(2), centre(2) + vecLen * curVect(2)];
        ZZ = [centre(3), centre(3) + vecLen * curVect(3)];

        plot3(XX, YY, ZZ, 'Color', curColor);
        hold on

    end
end





%                     leftInds  = sub2ind(size(marksPrev), indsYLeft,  indsXLeft);
%                     rightInds = sub2ind(size(marksPrev), indsYRight, indsXRight);
%                     lefts       = marksPrev(leftInds);
%                     rights      = marksPrev(rightInds);
%                     leftsMask = mask(leftInds);
%                     rightsMask = mask((rightInds));
%                     
%                     % check is something is empty
%                     lenEmpLeft = length(leftsMask(leftsMask == 0))/length(indsYLeft);
%                     lenEmpRight = length(rightsMask(rightsMask == 0))/length(indsYRight);
%                     
%                     if lenEmpLeft >= 0.5 % left should be an empty cell
%                         left = emptyCellID;
%                         depthLeft = depthCentral;
%                     else
%                         lefts = lefts(lefts>0);
%                         if isempty(lefts)
%                             left = 0;
%                         else
%                             left = mode(lefts);
%                             coordsLeft = [min(rows(j) + displacements(2,1), r), min(cols(j) + displacements(2,2),c)];
%                             depthLeft =  I(coordsLeft(1), coordsLeft(2));
%                         end
% 
%                     end
%                     
%                     if lenEmpRight >= 0.5 % left should be an empty cell
%                         right = emptyCellID;
%                         depthRight = depthCentral;
%                     else
%                         rights = rights(rights>0);
%                         if isempty(rights)
%                             right = 0;
%                         else
%                             right = mode(rights);
%                             coordRight = [min(rows(j) + displacements(3,1), r), min(cols(j) + displacements(3,2), c)];
%                             depthRight = I(coordRight(1), coordRight(2));   % learn from exact positions
%                         end  
%                     end














