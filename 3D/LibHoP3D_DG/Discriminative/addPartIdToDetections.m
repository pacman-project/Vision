load('statistics_3_2_7', '-mat');
load('statisticsAggregated_3_2_7','-mat');
output = [outputStatistics(:,2),outputStatistics(:,1),outputStatistics(:,4)];
partIds = zeros(length(output),1);
for i = 1:length(X)
    disp(i);
    inds = find(output(:,1) == X(i,1) & output(:,2) == X(i,2) & output(:,3) == X(i,3));
   for j=1:length(inds)
      partIds(inds(j)) = i; 
   end
   clear inds;
end