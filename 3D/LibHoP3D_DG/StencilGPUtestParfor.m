% gpuMode = 0   only CPU is used
% gpuMode = 1   only gpuArray is used
% gpuMode = 2   arrayfun + nested function
% gpuMode = 0   

function StencilGPUtestParfor

    function X = updateParentGrid(row, col, N)
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


    gpuMode = 2;

    gpu = gpuDevice(1);
    gridSize = 500;
    numGenerations = 100;
    initialGrid = (rand(gridSize,gridSize) > .75);


    % Draw the initial grid
    hold off
    imagesc(initialGrid);
    colormap([1 1 1;0 0.5 0]);
    title('Initial Grid');

    if gpuMode == 0
        grid = initialGrid;
    elseif gpuMode >0
        grid = gpuArray(initialGrid);
    end



    if gpuMode < 2
        timer = tic();
        % Loop through each generation updating the grid and displaying it
        for generation = 1:numGenerations
            grid = updateGrid(grid, gridSize);
        end

        if gpuMode > 0
            grid = gather(grid);
            wait(gpu);
        end

        cpuTime = toc(timer);
        fprintf('Average time : %2.3fms per generation.\n', ...
            1000*cpuTime/numGenerations);

        imagesc(gather(grid));
        title(num2str(generation));
        drawnow;


    elseif gpuMode == 2
        
        timer = tic();

        rows = gpuArray.colon(1, gridSize)';
        cols = gpuArray.colon(1, gridSize);
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

    end

end