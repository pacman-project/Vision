% this is to try implementing a function imdilate with GPU

function  ImDilateGPUTrial()

    function a1 = imdilateMy(rows, cols)
        sum = 0;
        for ii = rows - dy:rows+dy
            if ii > 0
                for jj = cols - dx:cols+dx
                    if jj > 0
                        if a(ii, jj) > 0
                            sum = sum + 1;
                        end
                    end
                end
            end
        end
        if sum
            a1(rows, cols) = 1; 
        end
    end

    pointerFunc = @imdilateMy;

    gpu = gpuDevice(1);
    is_GPU_USED = true;
    
    tic
        dx = 4;
        dy = 1;
        sizeA = 500;
        
        nm = [3,9];
        st = strel('rectangle', nm);
        for i = 1:10
            
            
            if is_GPU_USED
                a = gpuArray(uint8(round(rand(sizeA)*0.51)));
                a1 = gpuArray(uint8(zeros(sizeA, sizeA)));
            else
                a = uint8(round(rand(500)*0.51));
            end
            
%             a = imdilate(a, st);
%             str = ['Images/',num2str(i), '.png'];
%             a = a * 255;
%             if is_GPU_USED
%                 a = gather(a);
%             end
            [rowsC, colsC] = find(a > 0);
            
            arrayfun(pointerFunc, rowsC, colsC);
            
            a = a1; 
            imwrite(a, str, 'png');
            
        end
    
        wait(gpu);
        reset(gpu);
    toc
    
    



end

