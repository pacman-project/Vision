% this function detects the second layer parts, however it does not perform
% inhibition


function [nearestClusters, errors] = discretizeY(yMarks, yErrors, yMarks_alt, yErrors_alt, strLen, offsetIn) % prepare data for line discretization
    
        
    nearestClusters = zeros(strLen,1);
    errors = zeros(strLen,1);
    
    if strLen < 5
        return;
    end
    
    middle = 3;
    offset = 2;
    
    if offsetIn == 1
        start = 3;
        disp('ERROR 51');
        return;
    elseif offsetIn == 0  % should always be the case!
        start = 3;
    end
    
    
    for i = start:2:strLen-1
        central = yMarks(i, middle);
        
        tops =    yMarks(i-offset, middle - offset :2: middle + offset);
        bottoms = yMarks(i+offset, middle - offset :2: middle + offset);
        
        indsT = find(tops == central);
        lenT = length(indsT);
        indsB = find(bottoms == central);
        lenB = length(indsB);
        
        if lenT > 0 && lenB> 0
            nearestClusters(i) = central;
            errors(i) = (yErrors(i, middle) + yErrors(i-offset, indsT(1)) + yErrors(i+offset, indsB(1))) / 3;
        end

    end



end




%     for i = start:2:strLen-1
%         
%         line = yMarks(i-2:2:i+2, middle);
%         lineErr = yErrors(i-2:2:i+2, middle);
%         lineAlt = yMarks_alt(i-2:2:i+2, middle);
%         lineErrAlt = yErrors_alt(i-2:2:i+2, middle);
% 
% 
%         candidateEl = mode(line);     
%         inds = find(line==candidateEl);
%         kolEl = length(inds);  % the most frequent element
% 
% 
%         if kolEl == 1 
% 
%             if line(2) ~= 0
%                 candidateEl = line(2);
% 
%                 lineTop = yMarks(i-2, middle-offset:middle+offset);
%                 lineBottom = yMarks(i+2, middle-offset:middle+offset);
% 
%                 indsTop = find(lineTop==candidateEl);
%                 kolElTop = length(indsTop); 
%                 indsBottom  = find(lineBottom==candidateEl);
%                 kolElBottom = length(indsBottom);
% 
%                 if kolElBottom>0 && kolElTop>0
%                     nearestClusters(i) = candidateEl;
%                     errors(i) = (lineErr(inds(1)) + yErrors(i-2, indsTop(1)) + yErrors(i+2, indsBottom(1))) / 3;
%                 end
%             end
% 
% 
%         elseif kolEl == 2
%             nearestClusters(i) = candidateEl;
%             % check whether there is the same element in alternative marks
%             indsAlt = find(lineAlt == candidateEl);
%             kolAlt = length(indsAlt);
%             if kolAlt >= 1 % element is approved by the alternative hypothesis
%                 sumErrAlt = sum(lineErrAlt(indsAlt))/kolAlt; 
%                 errors(i) = (lineErr(inds(1)) + lineErr(inds(2)) + sumErrAlt) / 3;
%             else % not approved (error increases)
%                 errors(i) = 1.3 * (lineErr(inds(1)) + lineErr(inds(2))) / 2; % TODO better
%             end
%         elseif kolEl == 3;
%             nearestClusters(i) = candidateEl;
%             errors(i) = (lineErr(inds(1)) + lineErr(inds(2)) + lineErr(inds(3))) / 3;
%         end
%     end




