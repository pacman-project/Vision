% this is to analize layer 4 aggregated

stats = load('statistics/Layer4Aggregated.mat');
stats = stats.triples4Layer;

X = zeros(1000000,4);
ind = 0;

Z = load('statistics/layer3_final.mat');
Z = Z.adopted;  % these are 3rd layer elements

[rZ, cZ] = size(Z);


% rewrite stats as a feature vector

for i = 1:rZ
    for j = 1:rZ
        for k = 1:rZ
            freq = stats(i,j,k);
            if freq > 0
                ind = ind + 1;
                X(ind,:) = [i,j,k, freq];
            end
        end
    end
end

X = X(1:ind,:);