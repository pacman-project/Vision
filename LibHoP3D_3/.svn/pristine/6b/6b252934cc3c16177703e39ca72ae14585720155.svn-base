% this is cross-validation for SVM with gaussian kernels

function [C, accuracy] = RadialCrossValidation_Chi(X, y, fold, threads)

    if nargin<3
        fold = 5;
    end
    max = 0;
    bestC = 0;
    
    Cs = [0.5, 1, 5, 10, 50, 100, 500, 1000, 5000, 10000, 50000];
    
    lenC = length(Cs);  
    AccT = zeros(lenC, 1);
    
    for i = 1:lenC 

        curC = Cs(i);
        string = ['-t 5 -c ',num2str(curC), ' -v ', num2str(fold), ' -z ', threads];
        acc = svmtrain(y, X, string);
        AccT(i) = acc;
        if (acc > max)
            max = acc;
            bestC = curC;
        end

    end
    C = bestC;
    accuracy = max;
    AccT
end