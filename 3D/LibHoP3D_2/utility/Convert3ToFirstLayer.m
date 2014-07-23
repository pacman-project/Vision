% returns a table of size (lenCombs x lenX1)

function [XX] = Convert3ToFirstLayer(X, lenCombs, table2)

    %  convert X to 6 dimensional vector XX
    XX = zeros(lenCombs, 6);

    for i = 1:lenCombs
        curEl = X(i,:);

        left   = curEl(1);
        centre = curEl(2);
        right  = curEl(3);

        clusterXL = table2(left, 1);
        clusterYL = table2(left, 2);
        clusterXR = table2(right, 1);
        clusterYR = table2(right, 2);
        clusterXC = table2(centre, 1);
        clusterYC = table2(centre, 2); 

        XX(i, :) = [clusterXL, clusterYL, clusterXC, clusterYC, clusterXR, clusterYR];
    end
    
end