% this is to prepare secon layer slices to be fed to lineDiscretization
% function


function [nearestClusters, errors] = discretizeY(yMarks, yErrors, yMarks_alt, yErrors_alt, strLen, offset) % prepare data for line discretization

    nearestClusters = zeros(strLen,1);
    errors = zeros(strLen,1);
    
    if offset == 1
        start = 3;
    elseif offset == 0  % should always be the case!
        start = 3;
    end
    
    for i = start:2:strLen-1
        
        line = yMarks(i-2:2:i+2);
        lineErr = yErrors(i-2:2:i+2);
        lineAlt = yMarks_alt(i-2:2:i+2);
        lineErrAlt = yErrors_alt(i-2:2:i+2);
        
        candidateEl = mode(line);
        inds = find(line==candidateEl);
        kolEl = length(inds);
        
        if kolEl == 2
            nearestClusters(i) = candidateEl;
            % check whether there is the same element in alternative marks
            indsAlt = find(lineAlt == candidateEl);
            kolAlt = length(indsAlt);
            if kolAlt >= 1 % element is approved by the alternative hypothesis
                sumErrAlt = sum(lineErrAlt(indsAlt))/kolAlt; 
                errors(i) = (lineErr(inds(1)) + lineErr(inds(2)) + sumErrAlt) / 3;
            else % not approved (error increases)
                errors(i) = 1.5 * (lineErr(inds(1)) + lineErr(inds(2))) / 2;
            end
        elseif kolEl == 3;
            nearestClusters(i) = candidateEl;
            errors(i) = (lineErr(inds(1)) + lineErr(inds(2)) + lineErr(inds(3))) / 3;
        end
    end
end