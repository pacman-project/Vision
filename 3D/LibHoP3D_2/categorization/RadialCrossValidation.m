% this is cross-validation for SVM with gaussian kernels

function [C, G, accuracy] = RadialCrossValidation(X, y, fold, threads)
    if nargin<3
        fold = 5;
    end
    max = 0;
    bestC = 0;
    
    Cs = [0.1, 1, 10, 100, 1000, 10000, 100000];
    Gs = [0.0002, 0.001, 0.01, 0.1, 1];
    
    lenC = length(Cs);
    lenG = length(Gs);
    
    AccT = zeros(lenC, lenG);
    
    for i = 1:lenC 
        for j = 1:lenG
            curC = Cs(i);
            curG = Gs(j);
            string = ['-t 2 -c ',num2str(curC), ' -g ', num2str(curG), ' -v ', num2str(fold)];  % , ' -z', threads
            acc = svmtrain(y, X, string);
            AccT(i,j) = acc;
            if (acc > max)
                max = acc;
                bestC = curC;
                bestG = curG;
            end
        end
    end
    C = bestC;
    G = bestG;
    accuracy = max;
    AccT
end