% This is an alternative function to compute quantiles

function [Y1, Y2, YM] = quantileMy(X, quant1, quant2, quantMiddle)

    if quant1 > 1 || quant2 > 1 || quantMiddle > 1
        disp('wrong quantiles!!');
    end
    
    sumX = sum(X);
    if sumX == 0;
        disp('Error. Empty sequence');
    end
    q1 = sumX * quant1;
    q2 = sumX * quant2;
    qM = sumX * quantMiddle;

    curSum = 0;
    i = 0;
    while curSum < q1
        i = i+1;
        curSum = curSum + X(i);
    end
    Y1 = i;
    
    while curSum < qM
        i = i+1;
        curSum = curSum + X(i);
    end
    YM = i;
    
    while curSum < q2
        i = i+1;
        curSum = curSum + X(i);
    end
    Y2 = i;
end

