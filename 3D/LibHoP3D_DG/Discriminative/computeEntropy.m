% numCategories = 51;
% normalParFreq = zeros(numCategories,1);
% logNormalParFreq = zeros(numCategories,1);
% subSumEntropy = zeros(numCategories,1);
% for i = 1:length(X)
%     disp(i);
%         normalParFreq = partFreq(i,:) / sum(partFreq(i,:));
%         logNormalParFreq = log(normalParFreq);
%         subSumEntropy = normalParFreq
%     inds = find(partIds == i);
%         for j=1:length(inds)
%             partFreq(i,categoryIds(inds(j))) = partFreq(i,categoryIds(inds(j)))+1;
%         end
% end

partEntropy = zeros(length(partFreq),1);
for i = 1:length(partFreq)
     disp(i);
     partEntropy(i) = theirEntropy(partFreq(i,:));
end
%         normalParFreq = partFreq(i,:) / sum(partFreq(i,:));
%         logNormalParFreq = log(normalParFreq);
%         subSumEntropy = normalParFreq
%     inds = find(partIds == i);
%         for j=1:length(inds)
%             partFreq(i,categoryIds(inds(j))) = partFreq(i,categoryIds(inds(j)))+1;
%         end
% end
