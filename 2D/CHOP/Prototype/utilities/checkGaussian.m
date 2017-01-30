function [result, val] = checkGaussian(x, alpha)
%     numberOfBins = 10;
%     x = (x - mean(x));
%     buffZone = 1/(2*numberOfBins);
%     histCenters = (-1/2+buffZone):2*buffZone:(1/2-buffZone);
%     xHist = hist(x, histCenters);
% %    peakfinder(xHist, (max(xHist)-min(xHist))/numberOfBins, [], 1, true, true);
%     peaks = peakfinder(xHist, limit, [], 1, true, true);
%     result = numel(peaks) == 1;
    
    % Additional check. 
    val = std(x);
    if val > alpha
         result = false;
    else
         result = true;
    end
end