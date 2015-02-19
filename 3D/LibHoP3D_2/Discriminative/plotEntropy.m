function [  ] = plotEntropy( table, id )

    x = 1:51;
    y = table(id, :);
    
    plot(x, y);

end

