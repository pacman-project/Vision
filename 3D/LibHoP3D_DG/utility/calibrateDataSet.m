% this function uses images of a sphere or a cylinder to establish ratio
% between x and z dimensions

function zScale = calibrateDataSet(dataSetNumber, depthPath)
    
    if dataSetNumber == 3
        str = [depthPath, '/99_4_1_1_2_2.png'];
        I = imread(str);
        I = I(:,:,1);
        I = double(I);

        % determine radius of the sphere
        [r,c]= find(I>0);
        rMin = min(r);
        rMax = max(r);
        cMin = min(c);
        cMax = max(c);

        diam = rMax - rMin;
        rad = diam/2;
        I= I(rMin:rMax, cMin:cMax);

        centerX = round(diam/2);
        centerY = round(diam/2);

        dDepth = I(centerY, centerX) - I(centerY, centerX + round(diam/4));
        zScale = rad*(1 - sqrt(3)/2)/dDepth;
    end
    

end

