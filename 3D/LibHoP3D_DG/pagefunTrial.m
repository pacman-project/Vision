function pagefunTrial()

    numObjects = 1000;
    roomDimensions = [50 50 5]; % Length * breadth * height in meters
    
    Map = randomTransforms(numObjects);
    Robot = randomTransforms(1);
    
    gMap.R = gpuArray(Map.R);
    gMap.T = gpuArray(Map.T);
    gRobot.R = gpuArray(Robot.R);
    gRobot.T = gpuArray(Robot.T);
    
    gpuPagefunTime = gputimeit(@()pagefunTransform(gRobot, gMap));
    fprintf(['It takes %3.4f seconds on the GPU using pagefun ',...
    'to execute %d transforms.\n'], gpuPagefunTime, numObjects);
    
    cpuTime = timeit(@()loopingTransform(Robot, Map));
    fprintf('It takes %3.4f seconds on the CPU to execute %d transforms.\n', ...
        cpuTime, numObjects);
        a = 2;


    function Tform = randomTransforms(N)
        Tform.T = zeros(3, N);
        Tform.R = zeros(3, 3, N);
        for i = 1:N
            Tform.T(:,i) = rand(3, 1) .* roomDimensions';
            Tform.R(:,:,i) = orth(rand(3, 3));
        end
    end

    function Rel = loopingTransform(Robot, Map)
        Rel.R = zeros(size(Map.R), 'like', Map.R); % Initialize memory
        Rel.T = zeros(size(Map.T), 'like', Map.T); % Initialize memory
        for i = 1:numObjects
            Rel.R(:,:,i) = Robot.R' * Map.R(:,:,i);
            Rel.T(:,i) = Robot.R' * (Map.T(:,i) - Robot.T);
        end
    end

    function Rel = pagefunTransform(Robot, Map)
        Rel.R = pagefun(@mtimes, Robot.R', Map.R);
        Rel.T = Robot.R' * bsxfun(@minus, Map.T, Robot.T);
    end



end