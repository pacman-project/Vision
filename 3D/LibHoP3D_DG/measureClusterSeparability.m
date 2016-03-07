function quantLIs = measureClusterSeparability(D)

    len = size(D, 1);
    quantLIs = zeros(1, len);
    
    for i = 1:len-1
        dd = D(i, :);
        dd = dd(dd>0);
        
        dd = sort(dd, 'ascend');
        lenDD = length(dd);
        within = zeros(1, lenDD);
        between = zeros(1, lenDD);
        
        for j = 1:lenDD
            within(j) = sum(dd(1:j))/j;
            between(j) = sum(dd(j+1:end))/(len-j-1);
        end
        
        x = 1:lenDD;
        plot(x, within./between);
%         hold on
%         plot(x, between);
        
        a = 2;
    end

end

