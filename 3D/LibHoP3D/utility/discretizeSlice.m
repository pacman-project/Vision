% this is to prepare secon layer slices to be fed to lineDiscretization
% function

function [nearestClusters, errors] = discretizeSlice(sliceMarks, sliceErrors, strLen) % prepare data for line discretization

    nearestClusters = zeros(strLen,1);
    errors = zeros(strLen,1);
    
    for i = 2:2:strLen-1
        
        line = [sliceMarks(i-1,3), sliceMarks(i,3), sliceMarks(i+1,3)];
        lineEr = [sliceErrors(i-1,3), sliceErrors(i,3), sliceErrors(i+1,3)];
        
        candidateEl = mode(line);
        [r,c] = find(line==candidateEl);
        kolEl = length(find(r));
        
        if kolEl == 2
            nearestClusters(i) = candidateEl;
            errors(i) = 1.5 * ((lineEr(r(1), c(1)) + (lineEr(r(2), c(2)))/2));
        elseif kolEl == 3;
            nearestClusters(i) = candidateEl;
            errors(i) = 1.5 * ((lineEr(r(1), c(1)) + (lineEr(r(2), c(2))) + (lineEr(r(3), c(3)))/3));
        end
        
        if kolEl >= 2 % in both cases
            window = sliceMarks(i-1:i+1,:);
            kol2 = length(window == candidateEl);
            deltaKol = kol2 - kolEl;
            if deltaKol > 0
                errors(i) = errors(i) - errors(i)*0.1*deltaKol;
            end
        end


    end


%         if (sliceMarks(i-1,3) == sliceMarks(i+1,3)) && (sliceMarks(i-1,3) ~= 0) % the simplest case
%             nearestClusters(i) = sliceMarks(i-1,3);
%             errors(i) = (sliceErrors(i-1,3) + sliceErrors(i+1,3))/2;
%         elseif (sliceMarks(i-1,3) ~= 0) || (sliceMarks(i+1,3) ~= 0) % at least one is not zero
%             nEls = [sliceMarks(i-1,3), sliceMarks(i+1,3), sliceMarks(i-1,1), sliceMarks(i-1,5), sliceMarks(i+1,1), sliceMarks(i+1,5)];
%             candidateEl = mode(nEls);
%             inds = find(nEls == candidateEl);
%             if length(inds) == 2
%                 nearestClusters(i) = candidateEl;
%                 errors(i) = 0.5;
%             elseif length(inds) >= 2
%                 nearestClusters(i) = candidateEl;
%                 errors(i) = 1.0;
%             end
%         elseif (sliceMarks(i-1,3) == sliceMarks(i+1,3)) && (sliceMarks(i-1,3) ~= 0) % both are equal to zero
%             % we do not want empty cells
%             nearestClusters(i) = 5;
%             errors(i) = 100.0;    
%         end