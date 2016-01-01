function [partsV, is_ok] = VisualizePart(layerID, partID, position, Norm, Xtemp, Ytemp, computeAxis, circleRad, vectLen, isVisualization)

    if nargin == 9
        isVisualization = true;
    end

    % make sure they are all row vectors
    if size(Norm, 1) ~= 1
        Norm = Norm';
    end
    if size(Xtemp, 1) ~= 1
        Xtemp = Xtemp';
    end
    if size(Ytemp, 1) ~= 1
        Ytemp = Ytemp';
    end
    
    % make sure they are all normalized
    Norm = Norm/norm(Norm);
    if ~isempty(Xtemp)
        Xtemp = Xtemp/norm(Xtemp);
    end
    if ~isempty(Ytemp)
        Ytemp = Ytemp/norm(Ytemp);
    end

    if computeAxis % this makes Xtemp and Ytemp axis to be aligned with global axis)
        [Xtemp, Ytemp] = ComputePartAxis(Norm);
    end
 
    partsV = PreparePartForVisualization_VI(partID, layerID, position, Norm, Xtemp, Ytemp, vectLen);
    
    if isVisualization
        is_ok = VisualizePart_VI(partsV, size(partsV, 1), circleRad);
    else
        is_ok = true;
    end
end

function [Xtemp, Ytemp] = ComputePartAxis(Norm)
    Xtemp = [1,0,0];
    Ytemp = cross(Norm, Xtemp);
    Xtemp = cross(Ytemp, Norm);
    Xtemp = Xtemp/norm(Xtemp);
    Ytemp = Ytemp/norm(Ytemp);
end