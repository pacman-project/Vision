% computes distance of two quaternions 
% according to ALG 5 of Kuffner's paper

function dist = QuanternionDistance(Q1, Q2)
    dist = sqrt(1 - norm(dot(Q1, Q2)));
end

