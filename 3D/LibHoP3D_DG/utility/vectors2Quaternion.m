function q = vectors2Quaternion(v1, v2)
    q = zeros(4,1);
    q(1:3) = cross(v1, v2);
    q(4) = sqrt((norm(v1)^ 2) * (norm(v2) ^ 2)) + dot(v1, v2);
    q = q'/norm(q);
end