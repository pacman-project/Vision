function gpuTest

    size = 7000;
    gpu = gpuDevice(1);
    
    a = rand(size);
    b = rand(size);

%     a = gpuArray.rand(size);
%     b = gpuArray.rand(size);
    
    tic
    
        c = a*b;
%         c1 = gather(c);
        
    toc
    
    reset(gpu)

end

