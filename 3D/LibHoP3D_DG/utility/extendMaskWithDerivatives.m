function mask = extendMaskWithDerivatives(mask, cluster1Bounds, Ix, Iy)

    mask(Ix < cluster1Bounds(1)) = 0;
    mask(Ix > cluster1Bounds(end)) = 0;
    mask(Iy < cluster1Bounds(1)) = 0;
    mask(Iy > cluster1Bounds(end)) = 0;   
end

