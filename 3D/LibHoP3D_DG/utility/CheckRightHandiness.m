function  DF = CheckRightHandiness(DF)

    lenF = size(DF, 1);
    
    numRightHanded = 0;
    numLeftHanded = 0;    
    
    for i = 1:lenF
        
        DFCur = reshape(DF(i, :), [3,3]);
        N = DFCur(:,1);
        X = DFCur(:,2);
        Y = DFCur(:,3);  % local x and y axis
        
        temp = -dot(cross(X, Y), N);
        if temp > 1.01 || temp < 0.99
            numLeftHanded = numLeftHanded + 1;
        elseif -temp > 1.01 || -temp < 0.99
            numRightHanded = numRightHanded +1;
            
            DF(i, 1:3) = -DF(i, 1:3);
        end
        
    end
    
    disp(numRightHanded);
    disp(numLeftHanded);

end

