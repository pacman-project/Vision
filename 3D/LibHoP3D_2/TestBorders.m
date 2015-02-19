% this function is to test functions getDisplacements and
% checkImageBoundaries.

layerID = 3;
centerX = 5;
centerY = 15;
displ3 = 6;
elementType = 2; 
maxRadius = 3;

[indsXLeft, indsYLeft, indsXRight, indsYRight] = getDisplacements(layerID, centerX, centerY, displ3, elementType, maxRadius);

r = 20;
c = 20;
a = zeros(r, c);
[indsXLeft, indsYLeft, indsXRight, indsYRight] = checkImageBoundaries(indsXLeft, indsYLeft, indsXRight, indsYRight, r, c);

leftInds  = sub2ind(size(a), indsYLeft,  indsXLeft);
rightInds = sub2ind(size(a), indsYRight, indsXRight);


a(leftInds) = 1;
a(rightInds) = 1;

imtool(a);

%  Everithing works well
