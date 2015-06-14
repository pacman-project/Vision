r = 100;
c = 100;

I = zeros(r,c);
size = [30,20];

% define frame of reference
X = rand(3,1);
X = X/norm(X);
tY = [0,1,0]';
Z = cross(X,tY);
Y= cross(X,Z);
V = [Z,X,Y];



vecLen = 10;
vectColors =eye(3);


tempX = V(:,2)/norm(V(:,2));
tempY = V(:,3)/norm(V(:,3));

vectX = size(1) * tempX;
vectY = size(2) * tempY;


centre = [50,50,50];
I(centre(2), centre(1)) = 1;

% vectX_Im = vectX(1:2);
% vectY_Im = vectY(1:2);


points = [centre - vectX' + vectY'; centre + vectX' + vectY';  ...
          centre + vectX' - vectY'; centre - vectX' - vectY'];
      points = round(points);

scatter3(points(:,1), points(:,2), points(:,3));
hold on

for i = 1:4
    I(points(i, 2), points(i,1)) = 1; 
end

imtool(I);

% plot the eigenvectors in 3D
for ii = 1:3
    curVect = V(:, ii);
    curColor = vectColors(ii, :);
    XX = [centre(1), centre(1) + vecLen * curVect(1)];
    YY = [centre(2), centre(2) + vecLen * curVect(2)];
    ZZ = [centre(3), centre(3) + vecLen * curVect(3)];

    plot3(XX, YY, ZZ, 'Color', curColor);
    hold on
end

bw = poly2mask(points(:,1), points(:,2), r, c);
I(bw == 1) = 10;
% bw2 = poly2mask(points(2:4,1), points(2:4,2), r, c);
imtool(bw);

