% this script tests quaternion distances

Q1 = [0,0,0,1];
maxA = 10;

x = 1:maxA;
y = zeros(1,maxA);

V1 = [cosd(0), 0, sind(0)]';

for alpha = 1:maxA
    V = [cosd(alpha), 0, sind(alpha)]';
    Q = computeQuaternion(V1, V);
    Qd = QuanternionDistance(Q1, Q);
    y(alpha) = Qd;
end

plot(x,y);