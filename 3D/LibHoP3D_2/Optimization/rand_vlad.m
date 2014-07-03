% generates randomly either 1 or -1

function [ out ] = rand_vlad( )
    out = round(rand);
    if out == 0
        out = -1;
    end
end

