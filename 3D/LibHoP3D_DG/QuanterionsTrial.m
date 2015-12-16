% here we do some experiments with distances in 6D space (translation plus rotations)
% the objective is to measure a distance from the origin P = [0,0,0]; N = [0,0,1]
% and anothe point that is defined as P1 and normal N1;


% http://www.mathworks.com/matlabcentral/fileexchange/1176-quaternion-toolbox

% functions available
DCM = dcmfromeuler(pi/6,0,0, 'zyx') % converts from Euler vectors to matrix form

q = dcm2q(DCM)  % converts direction cosine matrices into quaternions.
% DCM1 = q2dcm(q) % from quaternion to a matrix

% qq = qconj(q); % conjugate of the quaternion
[v, phi] = qdecomp(q);

q = vectors2Quaternion(v1, v2); % quaternion that is required to rotate 
                                % one vector to another one
                                
qvrot(q, v1); % rotates vector according to quaternion

% TODO quaternions to convert one frame of reference to another one

a = 2;