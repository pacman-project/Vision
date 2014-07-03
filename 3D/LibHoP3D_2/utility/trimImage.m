% this function trims the image by mask

function [outImage, outMask] = trimImage(inImage, mask)       
    mr = max(mask');
    mc = max(mask);
    indr = find(mr);
    indc = find(mc);
    indrmin = indr(1);
    indrmax = indr(length(indr));
    indcmin = indc(1);
    indcmax = indc(length(indc));
    outMask = mask(indrmin:indrmax, indcmin:indcmax);
    outImage = inImage(indrmin:indrmax, indcmin:indcmax);
end