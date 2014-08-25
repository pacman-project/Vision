len = 10;

combs = zeros(1^3, len);
ind = 0;

numComb = 2^len; % number of possible combinations
combStart = 2^(len-2);

for i = combStart:numComb
     str = dec2bin(i-1, len);
     % check whether this is a good combination
     k = strfind(str, '00'); % features for bad combinations
     if isempty(k)
         kk = strfind(str, '1111');
         if isempty(kk) % combination looks good       
            curComb = zeros(1,len);
            for j = 1:len
                curComb(j) = str2double(str(j));
            end
            ind = ind + 1;
            combs(ind, :) = curComb;
         end
     end
end  

combs = combs(1:ind,:);

