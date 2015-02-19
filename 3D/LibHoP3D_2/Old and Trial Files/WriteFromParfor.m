function [] = WriteFromParfor()


    % DOES not work !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    stat = zeros(10^2, 5);
    cur = 0;
    last = 0;
    
    pointerWriteLine = @writeLine;
    
    tic
    parfor i = 1:10
        
        for j = 1:10
            a = rand(5,1);
            k = (i-1)*10+j;
            pointerWriteLine(a, k);
        end
        
    end
    
    function writeLine(aa, kk)
        disp(kk);
        last = kk;
        stat(kk,:) = aa;
    end
    
    toc
    
    a = 2;

end

