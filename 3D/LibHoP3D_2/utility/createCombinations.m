function [combs] = createCombinations(len)

    numComb = 2^len; % number of possible combinations
    combStart = 2^(len-2);
    combs = [];

    for i = combStart:numComb
         str = dec2bin(i-1, len);
         % check whether this is a good combination
         k = strfind(str, '000'); % features for bad combinations
         if isempty(k)
             kk = strfind(str, '11');
             if isempty(kk) % combination looks more or less good
                 curComb = zeros(1,len);
                 for j = 1:len
                     curComb(j) = str2double(str(j));
                 end
                 combs = [combs; curComb];
             end
         end
    end
    
end