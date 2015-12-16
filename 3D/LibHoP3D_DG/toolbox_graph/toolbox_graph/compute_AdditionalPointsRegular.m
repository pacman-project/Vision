% This function computes some points withing the triangle

function [Vadd, Nadd, pointIDx, numPoints, areCentral, pointIndexing] = compute_AdditionalPointsRegular(V, F, N, gridStep)%, noPointsIDs)

        
    lenF = size(F, 2);
    halfStep = gridStep/3;
    vecLen = 2*gridStep;
    sides = [1,2; 1,3; 2,3];
    otherPoints = [3,2,1];
    numPoints = zeros(1, lenF);
    areCentral = {};
    
%     if nargin == 4;
%         noPointsIDs = zeros(1, lenF);
%     end

    areas = compute_Areas(V, F);
    
    parfor k = 1:lenF

        pointsBase =  [V(:, F(1,k))'; V(:, F(2,k))'; V(:, F(3,k))'];
        NormalsBase = [N(:, F(1,k))'; N(:, F(2,k))'; N(:, F(3,k))'];
        
%         if noPointsIDs(k)   % add only ONE point, which is a weighted sum of three vertices
%             Normals = NormalsBase(1,:) +NormalsBase(2,:) + NormalsBase(3,:);
%             Normals = Normals/norm(Normals);
%             pointsGlob = (pointsBase(1,:) + pointsBase(2,:) + pointsBase(3,:))/ 3;
%             
%             Vadd{k} = pointsGlob;
%             Nadd{k} = Normals;
%             numPoints(k) = 1;
%             areCentral{k} = 1;
%             continue;
%         end
        
        % 1 take the largest side
        sideLengths = zeros(3,1);
        
        sideLengths(1) = sqrt(sum((pointsBase(sides(1,1), :) - pointsBase(sides(1,2),:)).^2));
        sideLengths(2) = sqrt(sum((pointsBase(sides(2,1), :) - pointsBase(sides(2,2),:)).^2));
        sideLengths(3) = sqrt(sum((pointsBase(sides(3,1), :) - pointsBase(sides(3,2),:)).^2));
        
        idx = find(sideLengths == max(sideLengths));
        if length(idx) > 1
            idx = idx(1);
        end
        
%         scatter3(pointsBase(:,1), pointsBase(:,2), pointsBase(:,3), 'green', 'marker', '.');
%         axis equal;
%         hold on

        % 2 take the largest side;
        tempX = pointsBase(sides(idx,1), :) - pointsBase(sides(idx,2),:);    
%         plotVect(tempX, 1, [1,0,0], pointsBase(sides(idx,2),:));
        b = pointsBase(otherPoints(idx),:) - pointsBase(sides(idx,2), :);
%         plotVect(b, 1, [1,0,0], pointsBase(sides(idx,2),:));
        
        projBLen = (tempX * b')/norm(tempX); % length of the projection of the vector b to tempX
        projB = tempX * projBLen/ norm(tempX);
        originTemp = pointsBase(sides(idx,2), :) + projB;
        tempY = b - projB;
        
%         plotVect(tempY, 1, [1,0,0], originTemp);
%         plotVect(tempX, 1, [1,0,0], originTemp);
        
        % 3 find the transformation from the one coordinate system to the other
        
        is_ok = true;
        
        if norm(tempX)~= 0 && norm(tempY) ~= 0
            T = eye(4); T(1,4) = -originTemp(1); T(2,4) = -originTemp(2); T(3,4) = -originTemp(3);
            R = eye(4,4);
            R(1:3, 1) = tempX/norm(tempX); 
            R(1:3, 2) = tempY/norm(tempY);
            R(1:3, 3) = [1;1;1];
            if abs(rcond(R)) > 10^-4 
                M = inv(R)*T; % from global coordinates to coordinates of the triangle
                M1 = inv(M);  % from the coordinates of the triangle to global ones
            else
                is_ok = false;
            end
        else
            is_ok = false;
        end

        cur = 0;
        if is_ok
            % estimate ranges in the triangle-based coordinate system
            P1 = M*[pointsBase(1,:),1]';
            P2 = M*[pointsBase(2,:),1]';
            P3 = M*[pointsBase(3,:),1]';

            xs = [P1(1), P2(1), P3(1)];
            ys = [P1(2), P2(2), P3(2)];

            maxX = max(xs);
            maxY = max(ys);
            minX = min(xs);
            points = zeros(100, 4);
            is_central = zeros(100,1);

            % points = [minX, 0, 0, 1;   maxX, 0,0,1;   0, maxY, 0,1];

            slope = maxY/minX;
            curI = 1;

            for i = minX+halfStep:  gridStep:  0 - halfStep
                curJ = 0;
                isSomething = false;
                for j = 0+halfStep: gridStep:  (minX - i)*slope-halfStep
                    curJ = curJ + 1;
                    cur = cur + 1;
                    points(cur, :) = [i,j,0,1];

                    if mod(curI, 9) == 4 && mod(curJ, 9) == 4
                        is_central(cur) = 3;
                    elseif mod(curI, 3) == 2 && mod(curJ, 3) == 2
                        is_central(cur) = 2;
                    else
                        is_central(cur) = 1;
                    end
                    isSomething = true;
                end
                if isSomething
                    curI = curI + 1;
                end
            end

            slope = maxY/maxX;
            curI = 1;
            
            for i = 0 + halfStep  :gridStep:  maxX - halfStep
                curJ = 0;
                isSomething = false;
                for j = (maxX - i)*slope - halfStep :-gridStep: 0 + halfStep
                    curJ = curJ + 1;
                    cur = cur + 1;
                    points(cur, :) = [i,j,0,1];
                    
                    if mod(curI, 9) == 4 && mod(curJ, 9) == 4
                        is_central(cur) = 3;
                    elseif mod(curI, 3) == 2 && mod(curJ, 3) == 2
                        is_central(cur) = 2;
                    else
                        is_central(cur) = 1;
                    end
                    isSomething = true;
                end
                if isSomething
                    curI = curI + 1;
                end
            end
        end
        
        if cur == 0 || ~is_ok  % add only ONE point, which is a weighted sum of three vertices
            Normals = NormalsBase(1,:) +NormalsBase(2,:) + NormalsBase(3,:);
            Normals = Normals/norm(Normals);
            pointsGlob = (pointsBase(1,:) + pointsBase(2,:) + pointsBase(3,:))/ 3;
            
            Vadd{k} = pointsGlob;
            Nadd{k} = Normals;
            numPoints(k) = 1;
            areCentral{k} = 3;
            continue;
        end

        points = points(1:cur,:);
        is_central = is_central(1:cur);
        
        if length(is_central(is_central == 2))< 2 % it should be a central point
            points(cur+1, :) = (P1 + P2 + P3)/3;
            is_central(cur + 1) = 2;
        end
        if length(is_central(is_central == 3))< 2
            points(cur+1, :) = (P1 + P2 + P3)/3;
            is_central(cur + 1) = 3;
        end
        
        pointsGlob = M1*points';
        pointsGlob = pointsGlob(1:3,:)';
        
        lenPG = size(pointsGlob, 1);
        
        %% compute surface normals in these points (baricentric coordinates)
        Sbig = areas(k);
        S3 = TriangleAreas(repmat(pointsBase(1,:), [lenPG,1]), repmat(pointsBase(2,:),[lenPG,1]), pointsGlob);
        S1 = TriangleAreas(repmat(pointsBase(2,:), [lenPG,1]), repmat(pointsBase(3,:),[lenPG,1]), pointsGlob);
        S2 = TriangleAreas(repmat(pointsBase(1,:), [lenPG,1]), repmat(pointsBase(3,:),[lenPG,1]), pointsGlob);
        
        w1 = S1/Sbig;
        w2 = S3/Sbig;
        w3 = S3/Sbig;
        Normals = w1*NormalsBase(1,:) + w2*NormalsBase(2,:) + w3*NormalsBase(3,:);
        lengths = sqrt(sum(abs(Normals).^2,2));
        Normals = Normals./repmat(lengths, [1, 3]);
           
        Vadd{k} = pointsGlob;
        Nadd{k} = Normals;
        numPoints(k) = cur;
        areCentral{k} = is_central;
              
%         scatter3(pointsGlob(:,1), pointsGlob(:,2), pointsGlob(:,3), 'red', 'marker', '.');
%         hold on
%         
%         ids = is_central == 2;
%         scatter3(pointsGlob(ids,1), pointsGlob(ids,2), pointsGlob(ids,3), 'blue', 'marker', 'o');
%         hold on
%         
%         ids = is_central == 3;
%         scatter3(pointsGlob(ids,1), pointsGlob(ids,2), pointsGlob(ids,3), 'green', 'marker', 'x');
%         hold on
%         
%         axis equal;     
%         a = 2;
    end
    numPointsOverall = sum(numPoints);
    pointIndexing = zeros(2, numPointsOverall);
    
    cur = 1;
    for i = 1:lenF
        pointIndexing(1, cur:cur + numPoints(i)-1) = i;
        pointIndexing(2, cur:cur + numPoints(i)-1) = 1:numPoints(i);
        pointIDx{i} = cur:cur + numPoints(i)-1;
        cur = cur + numPoints(i);
    end
    
    pointIndexing = pointIndexing';  
end

function plotVect(V, vecLen, vectColors, point)
    for ii = 1:1
        curVect = V(ii, :);
        curColor = vectColors(ii, :);
        XX = [point(1), point(1) + vecLen * curVect(1)];
        YY = [point(2), point(2) + vecLen * curVect(2)];
        ZZ = [point(3), point(3) + vecLen * curVect(3)];

        plot3(XX, YY, ZZ, 'Color', curColor);
        hold on
    end
end

function areas = TriangleAreas(P1, P2, P3)  % given three points

    lenP = size(P1, 1);

    a = vectDist(P1, P2);
    b = vectDist(P1, P3);
    c = vectDist(P2, P3);

    p = (a + b + c) / 2; % half of the perimeter
    areas = sqrt((p-a).*(p-b).*(p-c).*p); 
end

function dist = vectDist(V1, V2)
    dist = sqrt(sum((V1-V2).^2, 2));
end


% alternative solutions (not very good)

    %% generate weights  1
    % a = rand(1, numPoints);
    % b = rand(1, numPoints);
    % w1 = min(a, b);
    % w2 = abs(a - b);
    % w3 = 1 - max(a, b);

    %% generate weights  2
    % step = 0.1;
    % w = [];
    % 
    % for i = 0:step:1
    %     W1 = i;
    %     for j = 0:step:1-W1;
    %         W2 = j;
    %         W3 = 1 - W1 - W2;
    %         temp = [W1, W2, W3];
    %         w = [w; temp];
    %     end
    % end

