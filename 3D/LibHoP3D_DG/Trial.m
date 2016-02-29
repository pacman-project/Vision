clear all
close all

%% Points and angles (angle in radians)
angle = pi/6;
mypoint = [3 sqrt(3)]


%% Calculate 'center' from extrema
m = 2;     % width/2  (I want to rotate around the center of the image)
n = 1'    % height/2

%% Tx matrices
first = [...
    1 0 -m;
    0 1 -n;
    0 0 1];

third = [...
    1 0 m;
    0 1 n;
    0 0 1];

second = [...
    cos(angle) -sin(angle) 0;
    sin(angle) cos(angle) 0;
    0 0 1];

R = third* second* first
T = eye(3); T(1,3) = 2;
S = [2.5, 0,0,; 0,0.5,0; 0,0,1];

%% Use homogenous coords
mp_hom = [mypoint 1];

%% Calculate (note because we premultiply,
rotated_hom = T*S*R * mp_hom'
rotated = rotated_hom(1:2)'


points = [rotated; mypoint; [m,n]];
scatter(points(:, 1), points(:, 2));

