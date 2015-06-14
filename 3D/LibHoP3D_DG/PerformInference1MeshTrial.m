function PerformInference1MeshTrial()
    
    [V, F, ~] = meshRead('D:\Input Data\Meshes\Aim@Shape_all\D00001.obj');
    
    vecLen = 0.1;
    vectColors = [1,0,0];
    % visualize the mesh and some of the normals
    
%     trisurf(F,V(:,1),V(:,2),V(:,3), 'FaceColor', [0.1, 0.1, 0.1], 'EdgeColor', [0.3, 0.3, 0.3], 'FaceAlpha', 1.0);
%     light('Position',[-1.0,-1.0,100.0],'Style','infinite');
%     axis equal;
%     lighting phong;
%     hold on
%      
    numPoints = 10;

    [V,F] = check_face_vertex(V,F);

    lenV = size(V, 2);
    lenF = size(F, 2);
    
    threshPlanar = 8; % degrees

    fCentres = sum(F,1)/3;

    fring = compute_face_ring(F);
    vring = compute_vertex_ring(F);
    [N,Nf] = compute_normal(V,F);

    % get additional points in the middle (centre of gravity) of each triangle
    [VCent, NCent] = compute_gravitycentre(V, F, N);
%     areas = compute_Areas(V, F);
%     q30 = quantile(areas, 0.2);
%     q90 = quantile(areas, 0.9);
    maxNumPoints = 15;
    
%     numPoints = computeNumPoints(areas, q30, q90, maxNumPoints);
    numPoints = int32(10 * ones(lenF,1));
    [VAll, NAll] = compute_AdditionalPoints(V, F, N, numPoints, maxNumPoints);
    
    radii = zeros(1, lenF);

%     fid = 8500;
%     ids = F(:,fid);
%     
%     plotFrame(N(:, ids(1)), vecLen, vectColors, V(:, ids(1)));
%     plotFrame(N(:, ids(2)), vecLen, vectColors, V(:, ids(2)));
%     plotFrame(N(:, ids(3)), vecLen, vectColors, V(:, ids(3)));

%     for fid = 8000:5:85000
%         plotFrame(Nf(:, fid), vecLen, vectColors, VCent(:, fid));
%     end

    nf = 10; % can be 1:4,10

    parfor i = 1:lenF
        
        Vcentre = VCent(:, i);
        Ncentre = NCent(:, i);
        points = zeros(3, maxNumPoints*nf);  % VAll + V + Vcent
        Normals = zeros(3, maxNumPoints*nf);  % VAll + V + Vcent
        
        % Create a list of faces
        facesId = i;
        
        while length(facesId) < nf
            if length(facesId) == 1
                facesId = [facesId, fring{facesId}];
            else
                facesIdNew = facesId;
                for j = 1:length(facesId)
                    facesIdNew = union(facesIdNew, fring{facesId(j)});
                end
                facesId = facesIdNew;
            end
        end
        
        % EXTRACT POINTS AND NORMALS THAT BELONG TO THESE FACES
        for j = 1:nf
            startPoint = 1 + maxNumPoints * (j-1);
            endPoint = maxNumPoints * j;
            points(:,  startPoint:endPoint) = squeeze(VAll(facesId(j), :, :));
            Normals(:, startPoint:endPoint) = squeeze(NAll(facesId(j), :, :));      
        end
        
        % COMPUTE DISTANCES FROM THE CENTRAL POINT
        dists = (points - repmat(Vcentre, [1, length(points)])).^2;
        dists = sqrt(sum(dists,1));
        [dists, ids] = sort(dists, 'ascend');
        points = points(:, ids);
        Normals = Normals(:, ids);
        
        rad = 0;
        counter = 0;
        for j = 1:size(points, 2)
            if abs(angleVec(Ncentre, Normals(:,j))) <= threshPlanar
                rad = dists(j);
            else
                counter = counter + 1;
                
                if counter == 5
                    break;
                end
            end
        end
     
        radii(i) = rad;     
        
        if mod(i, 100) == 0
            i
        end
        
    end
    
%     radiis = sort(radii, 'ascend');
%     x = 1:lenF;
%     plot(x, radiis);
    
    % TRY TO VISUALIZE THE RECONSTRUCTION 
    V = zeros(4 * length(radii(radii > 0)), 3);
    F = zeros(2 * length(radii(radii > 0)), 3);
    j = 1; 
    
    radii = perform_mesh_smoothing(face,vertex,radii,options);
    
    for i = 1:length(radii)
        
        if radii(i) == 0
            continue;
        end
        
        position = VCent(:, i);
        Norm = NCent(:, i);
        Norm = Norm/norm(Norm);
        
        if dot(Norm, [1;0;0]) > 10^-4
            xtemp = [1;0;0];
        else
            xtemp = [0;1;0];
        end
        
        Ytemp = cross(Norm, xtemp);
        Xtemp = cross(Norm, Ytemp);
        Ytemp = (radii(i)/2)*(Ytemp/norm(Ytemp));
        Xtemp = (radii(i)/2)*(Xtemp/norm(Xtemp));
        Norm = (radii(i)/2)*Norm;
        
        % establish 4 vertices and two riangles here
        V((j-1)*4 + 1, :) = (position + Xtemp + Ytemp)';  % vertex
        V((j-1)*4 + 2, :) = (position + Xtemp - Ytemp)';  % vertex  
        V((j-1)*4 + 3, :) = (position - Xtemp - Ytemp)';  % vertex  
        V((j-1)*4 + 4, :) = (position - Xtemp + Ytemp)';  % vertex  
        
        F((j-1)*2 + 1, :) = [(j-1)*4 + 1, (j-1)*4 + 2, (j-1)*4 + 3]; % faces
        F((j-1)*2 + 2, :) = [(j-1)*4 + 1, (j-1)*4 + 3, (j-1)*4 + 4];
        
        j = j+1;
        
        if mod(i, 100) == 0
            i
        end
    end
    
    
    trisurf(F,V(:,1),V(:,2),V(:,3), 'FaceColor', [0.1, 0.1, 0.1], 'EdgeColor', [0.3, 0.3, 0.3], 'FaceAlpha', 0.4);
    light('Position',[-1.0,-1.0,100.0],'Style','infinite');
    axis equal;
    lighting phong;
    hold on
    
    a = 2;
end

function ed = euclDist(D1, D2)
    ed = sqrt(sum((D1 - D2).^2));
end

function an = angleVec(V1, V2) % angle of two vectors
    an = acos(dot(V1, V2) / sqrt(dot(V1, V1)*dot(V2, V2)));
    an = an * 180 / pi;
end

function plotFrame(V, vecLen, vectColors, point)
    for ii = 1:1
        curVect = V(:, ii);
        curColor = vectColors(ii, :);
        XX = [point(1), point(1) + vecLen * curVect(1)];
        YY = [point(2), point(2) + vecLen * curVect(2)];
        ZZ = [point(3), point(3) + vecLen * curVect(3)];

        plot3(XX, YY, ZZ, 'Color', curColor);
        hold on
    end
end

function numPoints = computeNumPoints(areas, q30, q90, maxPoints)
    
    for i = 1:length(areas)
        if areas(i) < q30
            numPoints(i) = 0;
        elseif areas(i) > q90
            numPoints(i) = maxPoints;
        else
            numPoints(i) = round(maxPoints*(areas(i) - q30)/(q90 - q30));
        end
    end

end



%     for i = 1:lenF
% 
%         is_planar = true;
% 
%         % set of vertices that can be fit with a planar patch
%         Vid = zeros(1,40);
%         Vcur = zeros(3,40);
%         
%         % initially we take points from 
%         Vid(1:3) = F(:, i);  
%         Vcur(:,1) = V(:,F(1, i));
%         Vcur(:,2) = V(:,F(2, i));
%         Vcur(:,3) = V(:,F(3, i));
%         
%         curCentre = Vadd(:, i);
%         normCentre = Nf(:, i);
%         planarRad = euclDist(curCentre, Vcur(:,3));
%         
%         candidates = []; % candidates to a planar patch
%         lv = 3;
% 
%         while is_planar
%             
%             if isempty(candidates)
%                 for j = 1:lv
%                     candidates = union(candidates, vring{Vid(j)}); 
%                 end
%             end
%             % circumcenter
%             
%             % sort all candidates by distances
%             dists = (V(:, candidates) - repmat(curCentre, [1, length(candidates)])).^2;
%             dists = sqrt(sum(dists,1));
%             
%             [dists, ids] = sort(dists, 'ascend');
%             candidates = candidates(ids);
%             
%             for j = 1:length(candidates) % sort candidates by distance
%                 if ismember(candidates(j), Vid)
%                     continue;
%                 else
%                     if abs(angleVec(N(:, candidates(j)), normCentre)) < threshPlanar
%                         planarRad = dists(j);
%                         lv = lv + 1;
%                         Vid(lv) = candidates(j);
%                     else
%                         is_planar = false;
%                         break;
%                     end
%                 end
%             end
%             
%             
%         end
% 
%         a = 2;
% 
%     %     V1 = V(:, curF(1));
%     %     V2 = V(:, curF(2));
%     %     V3 = V(:, curF(3));





%         for j = 1:5:numPoints*nf
%             plotFrame(Normals(:, j), vecLen, vectColors, points(:, j));
%         end
%    


