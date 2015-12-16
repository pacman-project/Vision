function partsV = PreparePartForVisualization_VI(partID, layerID, position, Norm, Xtemp, Ytemp, vectLen)
    
    if layerID == 2
        partsV = [position,Norm*vectLen];
        return;
    end
    
    positionCanonical = [0,0,0];
    OrientationCanonical = [0,0,1];

    [M,R,T] = TransformationMatrix(Xtemp, Ytemp, Norm, position);
    
    if layerID >= 3
        dd = load('Temp/Layer3/partSelection3.mat');
        lenOut{3} = dd.nNClusters{3};
        partsOut{3} = dd.partsOut;
        pairsAll{3} = dd.pairsAll;
    end
    if layerID >= 4
        dd = load('Temp/Layer4/partSelection4.mat');
        lenOut{4} = dd.nNClusters{4};
        partsOut{4} = dd.partsOut;
        pairsAll{4} = dd.pairsAll;
    end
    if layerID >= 5
        dd = load('Temp/Layer5/partSelection5.mat');
        lenOut{5} = dd.nNClusters{5};
        partsOut{5} = dd.partsOut;
        pairsAll{5} = dd.pairsAll;
    end
    if layerID >= 6
        dd = load('Temp/partSelection6.mat');
        lenOut{6} = dd.nNClusters{layerID};
        partsOut{6} = dd.partsOut;
        pairsAll{6} = dd.pairsAll;
    end
    
    if layerID == 3
        partsV = [pairsAll{3}(partsOut{3}(partID, 1), 3:end); [positionCanonical, OrientationCanonical]; pairsAll{3}(partsOut{3}(partID,3), 3:end)];
        for i = 1:size(partsV, 1)
            % convert positions of subparts
            point = [partsV(i, 1:3), 1];
            PointOut = M*point';
            partsV(i, 1:3) = PointOut(1:3)';

            % convert orientations of sub-parts
            orient = [partsV(i, 4:6), 1];
            orientOut = R*orient';
            partsV(i, 4:6) = vectLen * orientOut(1:3)/norm(partsV(i, 4:6));
        end
    end
    
    if layerID > 3
        partsV = [];
        curPart = partsOut{layerID}(partID, :);
        curPairs = pairsAll{layerID};
        % represent parts in CAnonical frame of reference
        partsV = PreparePartForVisualization_VI(curPart(2), layerID-1, positionCanonical, OrientationCanonical, [1,0,0], [0,1,0], vectLen);
        for i = 1:2:3 % two other subparts
            partPrevID = curPart(i);
            [Xt, Yt] = ComputePartAxis(curPairs(partPrevID, 6:8));
            partsV = [partsV; PreparePartForVisualization_VI(curPairs(partPrevID, 2), layerID-1, curPairs(partPrevID, 3:5), curPairs(partPrevID, 6:8), Xt, Yt, vectLen)];
        end
        
        % rotate and translate to the corre
        for i = 1:size(partsV, 1)
            % convert positions of subparts
            point = [partsV(i, 1:3), 1];
            PointOut = M*point';
            partsV(i, 1:3) = PointOut(1:3)';

            % convert orientations of sub-parts
            orient = [partsV(i, 4:6), 1];
            orientOut = R*orient';
            partsV(i, 4:6) = vectLen * orientOut(1:3)/norm(partsV(i, 4:6));
        end
    end
end

function [M,R,T] = TransformationMatrix(Xtemp, Ytemp, NormalCentral, central_pos)
    T = eye(4); T(1:3,4) = central_pos';
    R = eye(4,4); R(1:3, 1) = Xtemp; R(1:3, 2) = Ytemp; R(1:3, 3) = NormalCentral;
    M = T*R;
end

function [Xtemp, Ytemp] = ComputePartAxis(Norm)
    Xtemp = [1,0,0];
    Ytemp = cross(Norm, Xtemp);
    Xtemp = cross(Ytemp, Norm);
    Xtemp = Xtemp/norm(Xtemp);
    Ytemp = Ytemp/norm(Ytemp);
end


%     %% RotMatr(a, b) defines a transformation that rotates UNIT length vector a to b
%     GG = @(A,B) [ dot(A,B) -norm(cross(A,B)) 0;  norm(cross(A,B)) dot(A,B)  0;  0 0 1]; % a 2D roatational matrix
%     FFi = @(A,B) [ A (B-dot(A,B)*A)/norm(B-dot(A,B)*A) cross(B,A) ];
%     UU = @(Fi,G) Fi*G*inv(Fi);
%     U = UU(FFi(a,b), GG(a,b));
%     RotMatr = @(a,b) UU(FFi(a,b), GG(a,b));

