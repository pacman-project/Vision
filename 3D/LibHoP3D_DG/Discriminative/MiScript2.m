% load('TableModelLayer3.mat');
load('TableLayer3.mat');
table = table';
numCombs = size(table, 1);

denominator = sum(table, 1);
multiplier = max(denominator);
denom = repmat(denominator, [numCombs, 1]);

table = round(multiplier*(table ./ denom));



partEntropy = zeros(numCombs,1);

parfor i = 1:numCombs
     partEntropy(i) = theirEntropy(table(i, :));
     
     if mod(i, 10) == 0
         i
     end
end




