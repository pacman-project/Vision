% parce princeton Partitioning in training/test
empty = 9999;

curY = '';

Ys = [];
Xs = [];
cur = 0;

for i = 4:length(PSB1);
    
    if PSB1(i) == empty && PSB1(i+1) ~= empty
        curY = PSB{i};    
    elseif PSB1(i) ~= empty
        cur = cur+1;
        Ys{cur} = curY;
        Xs = [Xs; PSB1(i)];
    end
    
    
end


