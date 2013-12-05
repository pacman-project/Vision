% this is to display layer 4 elements one after another

Z = load('Statistics/layer4.mat');
Z = Z.Z;

Z = round(Z);
[r,c] = size(Z);
Z(Z<2) = 2;
Z(Z>8) = 8;

% pre-compute a table
table = zeros(9,9);

for i = 1:9
    for j = 1:9
        table(i,j) = (i - 1)* 9 + j;
    end
end

for i = 40:70
    curLine = zeros(1,9);
    for j = 1:9
        curLine(j) = table(Z(i, 2*j-1), Z(i,2*j));
    end
    figure;
    VisualizerLayer4(curLine);
end
        

