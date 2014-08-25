% Here we aggregate the first layer statistics, i.e. compute frequencies of
% all combinations

clear all;

nClusters = 9;
n2Clusters = nClusters * nClusters;

% cluster3stats = load('cluster3stats.mat');
% cluster3stats = cluster3stats.cluster3stats;

%triple3stats = load('statistics/Statistics5D_4layer_19_2300.mat');
%triple3stats = triple3stats.tripleStatistics;

triple3stats = load('statistics/Cluster_196.mat');
triple3stats = triple3stats.stats1;

% transform the file with fourth layer elements

[r,c] = size(triple3stats);
ind = 0;
subsetBegin = 1;

% c has to be equal to 9
tripleNew = zeros(100000,c+1);
lenNew = 0;

% first sort the matrix triple3stats by the first column
a = triple3stats(:,1);
[temp, inds] = sort(a);
triple3stats = triple3stats(inds, :);

for ii = 1:81
    
    % extract the subset
    subsetLen = length(a(a==ii));
    subset = triple3stats(subsetBegin:subsetBegin + subsetLen - 1, :);
    % sort by the second column to speed up
    subsetBegin = subsetBegin + subsetLen; % for the next subset (if any)
    
    % sort the extracted subset by the second column
    aa = subset(:,2);
    [temp, inds] = sort(aa);
    subset = subset(inds, :);
   
    % go through the subset and aggregate statistics
    for i = 1:subsetLen

        % read the current line
        line = subset(i,:);
        secondEl = line(2);
        if secondEl == 0
            continue;
        end
        % find number of the same lines in the data
        score = 1;

        for j = i+1:subsetLen   % count how many same lines exist in the data
            curSecondEl = subset(j,2);
            if curSecondEl == 0
                continue;
            elseif curSecondEl > secondEl
                break;
            elseif curSecondEl == secondEl

                nextLine = subset(j,:);
                if nextLine == line
                    score = score + 1;
                    subset(j,:) = zeros(1, c);
                end
            end
        end

        ind = ind + 1;
        tripleNew(ind, 1:c) = line;
        tripleNew(ind, c+1) = score;

        if mod(i,10) == 0
            i
        end

    end
end

tripleNew = tripleNew(1:ind,:);
%save('statistics/Statistics5D_4layer_19_2300_aggregated.mat','tripleNew');
save('statistics/Cluster_196_aggregated.mat','tripleNew');


% O(n^2) algorithm

% for i = 1:r
%     % read the current line
%     line = triple3stats(i,:);
%     if line(1) == 0
%         continue;
%     end
%     % find number of the same lines in the data
%     score = 1;
%     for j = i+1:r        
%         if triple3stats(j,1) == 0
%             continue;
%         end
%         nextLine = triple3stats(j,:);
%         
%         if nextLine == line
%             score = score + 1;
%             triple3stats(j,:) = [0,0,0,0,0];
%         end
%     end
%     ind = ind + 1;
%     tripleNew(ind, 1:5) = line;
%     tripleNew(ind,6) = score;
%     triple3stats(i,:) = [0,0,0,0,0];
%     
%     
%   %  if mod(i,10) == 0
%         i
%  %   end
% end
