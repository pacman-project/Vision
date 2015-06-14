% gpuMode = 0   only CPU is used
% gpuMode = 1   only gpuArray is used
% gpuMode = 2   arrayfun + nested function
% gpuMode = 0   

function ArrayFuncDistTrial1

    function X = updateParentGrid(row, col, N)
        P = PP(row, col);
        
        % Take account of boundary effects
        rowU = max(1,row-1);  rowD = min(N,row+1);
        colL = max(1,col-1);  colR = min(N,col+1);
        % Count neighbors
        neighbors ...
            = grid(rowU,colL) + grid(row,colL) + grid(rowD,colL) ...
            + grid(rowU,col)                   + grid(rowD,col) ...
            + grid(rowU,colR) + grid(row,colR) + grid(rowD,colR);
        % A live cell with two live neighbors, or any cell with
        % three live neighbors, is alive at the next step.
        X = (grid(row,col) & (neighbors == 2)) | (neighbors == 3);
    end


    gpu = gpuDevice(1);
    gridSize = 500;
    numGenerations = 100;
    initialGrid = (rand(gridSize,gridSize) > .75);


    % Draw the initial grid
    hold off
    imagesc(initialGrid);
    colormap([1 1 1;0 0.5 0]);
    title('Initial Grid');
    

    grid = gpuArray(initialGrid);

    timer = tic();

    rows = gpuArray.colon(1, gridSize)';
    cols = gpuArray.colon(1, gridSize);
    
    
    
    dim = 1000;
    A = gpuArray.rand(dim);
    
    lenP = 500;
    lenQ = 500;
    
    PP = gpuArray.rand(lenP, dim);
    QQ = gpuArray.rand(lenQ, dim);
    
    
    for generation = 1:numGenerations
        grid = arrayfun(@updateParentGrid, rows, cols, gridSize);
    end

    wait(gpu); % Only needed to ensure accurate timing
    gpuArrayfunTime = toc(timer);

    % Print out the average computation time and check the result is unchanged.
    fprintf('Average time : %2.3fms per generation.\n', ...
        1000*gpuArrayfunTime/numGenerations);

    imagesc(gather(grid));
    title(num2str(generation));
    drawnow;
    
    reset(gpu);



end