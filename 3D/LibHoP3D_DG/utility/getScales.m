% this function returns a set of scales for image downsampling

function [scales, lineAdders] = getScales(dataSetNumber, is_multyScale)

    if ~is_multyScale % single scale only
        scales = 1.0;
        lineAdders{1} = '1.00';
        return;
    end

    if dataSetNumber == 3  % Vladislav_STD
        scales = [1.0, 0.82 0.66, 0.50, 0.33];
    end
    
    for i = 1: length(scales)
        
        line = num2str(scales(i), '%f');
        dotPos = strfind(line, '.');
        line = line(1:dotPos+2);
        lineAdders{i} = line;

    end
end

