function [partsV, circleRad] = PreparePartForVisualization_VI(partID, layerID, position, Norm, Xtemp, Ytemp, vectLen)
    
    positionCanonical = [0,0,0];
    NormCanonical = [0,0,1];
    frameCanonical = eye(3);
    OrientationCanonicalQ = computeQuaternion(NormCanonical, NormCanonical)';
    OrientationCanonicalQ246 = computeQuaternion(NormCanonical, Norm)';

    if layerID == 2
        
        partsV = [position, OrientationCanonicalQ246];  %  Norm*vectLen];
        circleRad = 1;
        return;
    elseif layerID == 4 && partID == 1
        partsV = [position, OrientationCanonicalQ246];  % Norm*vectLen];
        circleRad = 3;
        return
    elseif layerID == 4 && partID > 1
        partID = partID - 1;
    elseif layerID == 6 && partID == 1
        partsV = [position, OrientationCanonicalQ246];  % Norm*vectLen];
        circleRad = 9;
        return;
    elseif layerID == 6 && partID > 1
        partID = partID - 1;
    end


    [M,R,T] = TransformationMatrix(Xtemp, Ytemp, Norm, position);
    
    global partsOut;
    global pairsAll;
    
    if layerID == 3
        partsV = [pairsAll{3}(partsOut{3}(partID, 1), 3:9); [positionCanonical, OrientationCanonicalQ]; pairsAll{3}(partsOut{3}(partID,3), 3:9)];
        for i = 1:size(partsV, 1)
            % convert positions of subparts
            point = [partsV(i, 1:3), 1];
            PointOut = M*point';
            partsV(i, 1:3) = PointOut(1:3)';

            % convert orientations of sub-parts
            orientQ = partsV(i, 4:7);
            QGlob = dcm2q([Xtemp; Ytemp; Norm]);
            partsV(i, 4:7) = qmult(QGlob, orientQ);

%             orient = partsV(i, 4:6);
%             orientOut = R*orient';
%             partsV(i, 4:6) = vectLen * orientOut(1:3)/norm(partsV(i, 4:7));
        end
        circleRad = ones(size(partsV, 1), 1);
    end
    
    if layerID >= 4
        curPart = partsOut{layerID}(partID, :);
        curPairs = pairsAll{layerID};
        % represent the central parts in canonical frame of reference
        
        [partsV, circleRad]           = PreparePartForVisualization_VI(curPart(2),              layerID-1, positionCanonical, NormCanonical,   [1,0,0], [0,1,0], vectLen);
        for i = 1:2:3 % two other subparts
            partPrevID = curPart(i);
            Q = curPairs(partPrevID, 6:9);
            Q = qnorm(Q);
            if layerID == 4
                NormC = qvrot(Q, NormCanonical);
                [Xt, Yt] = ComputePartAxis(NormC);
            else
                [Xt, Yt, NormC] = ComputeFrameRotated(frameCanonical, Q);
            end
            [partsVcur, circleRadCur] = PreparePartForVisualization_VI(curPairs(partPrevID, 2), layerID-1, curPairs(partPrevID, 3:5), NormC, Xt, Yt, vectLen);
            partsV = [partsV; partsVcur];
            circleRad = [circleRad; circleRadCur];
        end
        
        % rotate and translate
        for i = 1:size(partsV, 1)
            % convert positions of subparts
            point = [partsV(i, 1:3), 1];
            PointOut = M*point';
            partsV(i, 1:3) = PointOut(1:3)';

            % convert orientations of sub-parts
            orientQ = partsV(i, 4:7);
            QGlob = dcm2q([Xtemp; Ytemp; Norm]);
            partsV(i, 4:7) = qmult(QGlob, orientQ);
        end
    end
    
end

function [Xtemp, Ytemp, Norm] = ComputeFrameRotated(frameCanonical, Q)

    Xtemp = qvrot(Q, frameCanonical(1, :));
    Ytemp = qvrot(Q, frameCanonical(2, :));
    Norm =  qvrot(Q, frameCanonical(3, :));
end

function [M,R,T] = TransformationMatrix(Xtemp, Ytemp, NormalCentral, central_pos)
    T = eye(4); T(1:3,4) = central_pos';
    R = eye(4,4); R(1:3, 1) = Xtemp; R(1:3, 2) = Ytemp; R(1:3, 3) = NormalCentral;
    M = T*R;
end

function [Xtemp, Ytemp] = ComputePartAxis(Norm)
    Xtemp = [1,0,0];
    if abs(Xtemp*Norm') > 0.9
        Xtemp = [0,1,0];
        if abs(Xtemp*Norm') > 0.9
            Xtemp = [0,0,1];
        end
    end
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

