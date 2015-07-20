% this is to try geneartion of three random numbers
clear all;
pointsBase = [ 0.4636   -0.0643    0.0301;  0.4822   -0.0668    0.0316;  0.4644   -0.0602    0.0352];

gridStep = 0.001;
halfStep = gridStep/3;

% points = zeros(numPoints, 2);

%% generate weights  1
% a = rand(1, numPoints);
% b = rand(1, numPoints);
% w1 = min(a, b);
% w2 = abs(a - b);
% w3 = 1 - max(a, b);

%% generate weights  2
% step = 0.1;
% w = [];
% 
% for i = 0:step:1
%     W1 = i;
%     for j = 0:step:1-W1;
%         W2 = j;
%         W3 = 1 - W1 - W2;
%         temp = [W1, W2, W3];
%         w = [w; temp];
%     end
% end

%% generate weights  3

% 1 take the largest side
sides = [1,2; 1,3; 2,3];
otherPoints = [3,2,1];
sideLengths = zeros(3,1);
for i = 1:3
    sideLengths(i) = ComputeEuclDists(pointsBase(sides(i,1), :), pointsBase(sides(i,2),:), 3);
end
idx = find(sideLengths == max(sideLengths));
if length(idx) > 1
    idx = idx(1);
end

% 2 take the largest side;
tempX = pointsBase(sides(idx,1), :) - pointsBase(sides(idx,2),:);
lenBase = sqrt(sum(tempX.^2));
tempX = tempX/norm(tempX);
b = pointsBase(otherPoints(idx),:) - pointsBase(sides(idx,2), :);
projBLen = dot(tempX, b); % length of the projection of the vector b to tempX
projB = tempX * projBLen;
originTemp = pointsBase(sides(idx,2), :) + projB;
tempY = b - projB;
tempZ = [0,0,1];

% 3 find the transformation from the one coordinate system to the other

T = eye(4); T(1,4) = -originTemp(1); T(2,4) = -originTemp(2); T(3,4) = -originTemp(3);
R = eye(4,4); R(1:3, 1) = tempX; R(1:3, 2) = tempY/norm(tempY); R(1:3, 3) = [1;1;1];
M = inv(R)*T; % from global coordinates to coordinates of the triangle
M1 = inv(M);  % from the coordinates of the triangle to global ones

% estimate ranges in the triangle-based coordinate system
P1 = M*[pointsBase(1,:),1]';
P2 = M*[pointsBase(2,:),1]';
P3 = M*[pointsBase(3,:),1]';

xs = [P1(1), P2(1), P3(1)];
ys = [P1(2), P2(2), P3(2)];

maxX = max(xs);
maxY = max(ys);
minX = min(xs);
points = zeros(100, 4);
cur = 0;

% points = [minX, 0, 0, 1;   maxX, 0,0,1;   0, maxY, 0,1];

slope = maxY/minX;
for i = minX+halfStep:  gridStep:  0 - halfStep
    for j = 0+halfStep: gridStep:  (minX - i)*slope-halfStep
        cur = cur + 1;
        points(cur, :) = [i,j,0,1];
    end
end

slope = maxY/maxX;
for i = 0 + halfStep: gridStep: maxX - halfStep
    for j = (maxX - i)*slope - halfStep: -gridStep :0 + halfStep
        cur = cur + 1;
        points(cur, :) = [i,j,0,1];
    end
end

scatter(points(:,1), points(:,2), 'marker', '.');
axis equal;

pointsGlob = M1*points';
pointsGlob = pointsGlob';

% % compute points and visualize the plot
% numPoints = size(w,1);
% for i = 1:numPoints
%     points(i, :) = pointsBase(1, :) * w1(i) + pointsBase(2, :) * w2(i) + pointsBase(3, :) * w3(i);
% end

scatter3(pointsGlob(:,1), pointsGlob(:,2), pointsGlob(:,3), 'marker', '.');
a = 2;

