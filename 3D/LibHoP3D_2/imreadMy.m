function [ I, is_ok ] = imreadMy(str)

    is_ok = 1;
    try     
        I = imread(str);
    catch err
        try
            I = imread(str);
        catch err
            I = [];
            is_ok = 0;
        end
        disp('ERROR');
        disp(str);
    end


end

