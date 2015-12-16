% this function finds a quaternion which rotates vector V1 to V2

function Q = computeQuaternion(v1, v2)
    a = cross(v1, v2);
    w = norm(v1) * norm(v2) + v1'*v2;
    Q = qnorm([a; w]);
end

