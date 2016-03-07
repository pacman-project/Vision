function points = circleSampling(centre, rad, Norm, NumRadialSteps)

    if NumRadialSteps ~= 0
        adders = 1/NumRadialSteps;
    end

    if nargin == 4
        angleSteps = [45, 30, 20, 15, 10, 6, 5, 3, 3, 3, 2, 2, 2, 2, 2, 2];
    end

    tempX = [1,0,0];
    
    if size(Norm, 1) ~= 1  % Norm should be a row vector
        Norm = Norm';
    end
    if size(centre, 1) ~= 1  % Norm should be a row vector
        centre = centre';
    end
    
    if isequal(Norm'*tempX, Norm)  % normX == tempX 
        tempX = [0,1,0];
    end
    Y = cross(Norm, tempX);
    X = cross(Y, Norm);
    X = X/norm(X);
    Y = Y/norm(Y);
    
    numPoints = 1 + 2 * (round(360/angleSteps(end)));
    points = zeros(numPoints, 3);
    points(1, :) = centre;
    curPoint = 1;
    curAdder = 0;
    for j = 1:NumRadialSteps
        curAdder = curAdder + adders;
        curRad = rad *curAdder;
        for i = 2:(round(360/angleSteps(j))) + 1
            curPoint = curPoint + 1;
            points(curPoint, :) = centre + curRad * X * sind(angleSteps(j)*i) + curRad * Y * cosd(angleSteps(j)*i);
        end
    end
    points = points(1:curPoint, :);
    
    
end

