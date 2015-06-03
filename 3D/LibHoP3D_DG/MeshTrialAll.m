% This is to try everithing with meshes

[V, F, ~] = meshRead('D:\Input Data\Meshes\Aim@Shape_all\D00003.obj');
[V,F] = check_face_vertex(V,F);

lenV = size(V, 2);
lenF = size(F, 2);

% compute areas of all the triangles and angles
Ss = zeros(NumF, 1);

for i = 1:NumF
        
    curF = F(:,i);
    V1 = V(:, curF(1));
    V2 = V(:, curF(2));
    V3 = V(:, curF(3));

    a = V2 - V1;
    b = V3 - V1;
    c = V2 - V3;

    lenA = sqrt(dot(a,a));
    lenB = sqrt(dot(b,b));
    lenC = sqrt(dot(c,c));

    p = (lenA + lenB + lenC)/2;
    S = (p-lenA)*(p-lenB)*(p-lenC)*p;
    Ss(i) = S;

end

indMin = find(Ss == min(Ss));  % 2520
indMax = find(Ss == max(Ss));  % 11197

% Ss = sort(Ss, 'ascend');
% x = 1:lenF;
% plot(x, Ss);


a = 2;