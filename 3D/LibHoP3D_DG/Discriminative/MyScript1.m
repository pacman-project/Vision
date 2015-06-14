load('ready.mat') % 'partIDs', 'categoryIds', 'modelIds')
load('statisticsAggregated_3_2_7','-mat');

numParts = size(X,1);

lenStat = length(partIDs);

table = zeros(20, numParts);

for i = 1:lenStat
    
    if partIDs(i) > 0
        table(modelIds(i), partIDs(i)) = table(modelIds(i), partIDs(i)) + 1;
    end
    
    if mod(i, 100000) == 0
        i
    end
end

a = 2;