% this is Vladislav's GPU implementation of a function is_member

function is_memberMy()

    function isTrue = isMemberCur(cl1, cl2, cl3)
        
        isTrue = false;
        for j = 1:lenB
            isTrue = isTrue | (b(j,1) == cl1 & b(j,2) == cl2 & b(j,3) == cl3);
        end
    end

    len = 10^6;
    a = gpuArray.randi(100,len, 3);
    
    pointerisMemberCur = @isMemberCur;
    
    lenB = 20;
    
    for i = 1:12
        
        i
    
        b = randi(100,lenB, 3);

        tic

        ids = arrayfun(pointerisMemberCur, a(:,1), a(:,2), a(:,3));

        toc

        disp('My implementation:')
        str = ['number of lines matched: ', num2str(sum(ids))];
        disp(str);



        tic

        ids = ismember(a, b, 'rows');

        toc

        disp('Matlabs is_member function:')
        str = ['number of lines matched: ', num2str(sum(ids))];
        disp(str);
        
        disp('   ');
        
    end
    
    

    
    

end

