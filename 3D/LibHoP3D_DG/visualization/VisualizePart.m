% q is a quaternion

function [partsV, circleRads] = VisualizePart(layerID, partID, position, Q, Norm, Xtemp, Ytemp, computeAxis, circleRadBasic, vectLen, isVisualization)

    if nargin == 10
        isVisualization = true;
    end

    NormCanonical = [0,0,1];
    frameCanonical = eye(3);
    
    % make sure they are all row vectors
    if size(Q, 1) ~= 1
        Q = Q';
    end
    
    if ~computeAxis
        if size(Xtemp, 1) ~= 1
            Xtemp = Xtemp';
        end
        if size(Ytemp, 1) ~= 1
            Ytemp = Ytemp';
        end
        if size(Norm, 1) ~= 1
            Norm = Norm';
        end
        if ~isempty(Xtemp)
            Xtemp = Xtemp/norm(Xtemp);
        end
        if ~isempty(Ytemp)
            Ytemp = Ytemp/norm(Ytemp);
        end
        if ~isempty(Norm)
            Norm = Norm/norm(Norm);
        end
    end

    if layerID < 5
        if computeAxis
            % compute norm using quaternion
            Q = qnorm(Q);
            Norm = qvrot(Q, NormCanonical);
            [Xtemp, Ytemp] = ComputePartAxis(Norm);  % this makes Xtemp and Ytemp axis to be aligned with global axes
        end
    elseif layerID >= 5
        if computeAxis
            Q = qnorm(Q);
            [Xtemp, Ytemp, Norm] = ComputeFrameRotated(frameCanonical, Q);
        end
    end

    [partsV, circleRads] = PreparePartForVisualization_VI(partID, layerID, position, Norm, Xtemp, Ytemp, vectLen);
    
    if isVisualization
        is_ok = VisualizePart_VI(partsV, size(partsV, 1), circleRadBasic * circleRads, 0.01);
    else
        is_ok = true;
    end
end

function [Xtemp, Ytemp, Norm] = ComputeFrameRotated(frameCanonical, Q)

    Xtemp = qvrot(Q, frameCanonical(1, :));
    Ytemp = qvrot(Q, frameCanonical(2, :));
    Norm =  qvrot(Q, frameCanonical(3, :));
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