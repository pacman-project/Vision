% for this function we assume that centers of clusters are at equal
% distances from each other

function [output] = lineDiscretizationY(fx, strlen, kernelW, nClusters)

    % this is to compute score of an element to be at this position
    scores = zeros(1,strlen);
    output = zeros(1,strlen);
    
    for i = 3:strlen-2
        central = fx(i);
        window = fx(i-2:i+2);
        kol = length(find(window == central));
        kolSmaller = length(find(window == central-1));
        kolLarger = length(find(window == central+1));
        kolAdd = kolSmaller + kolLarger;
        scores(i) = kol + 0.5 * kolAdd;
        if scores(i) < 4
            scores(i) = 0;
        end
        
        % TODO todo something more smart with scores!
    end
    
    % performing discretization step here
    cur = 0;
    while (cur < strlen - 7)
        cur = cur + 3;
        % look for the max value in three rows (cur, cur+1, cur+2, ... , cur+5)
        matr =  scores(cur:cur+4);
        matr = matr .* kernelW;

        [c] = find(matr == max(matr));
        % if there are two similar maximums (not often happens!)
        if length(c) > 1
            c = c(length(c));
        end
        
        col = cur + c - 1;
        output(col) = fx(col);
        cur = col;        
    end 
    

   a = 10; 

end