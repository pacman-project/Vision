% this function finds a quaternion which rotates vector V1 to V2

function Q = computeQuaternion(v1, v2)

    if size(v1,1) == 1
        v1 = v1';
    end
    if size(v2,1) == 1
        v2 = v2';
    end
    a = cross(v1, v2);
    w = norm(v1) * norm(v2) + v1'*v2;

    Q = qnorm([a; w]);
end

