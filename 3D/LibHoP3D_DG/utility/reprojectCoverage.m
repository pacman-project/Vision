% this function reprojects coverage from the smaller to a higher scale
% this is done by replacing the 3rd channel of the hierarchy

function reprojectCoverage(list_input, lineAdder, lineAdderNext)
    
    lenF = length(list_input);
    
    for i = 1:lenF
        
        % read an image of the smaller resolution
        I = imread(list_input{i});
        I3 = I(:,:,3);
        
        % read the image of higher scale
        fileName = strrep(list_input{i}, lineAdder, lineAdderNext);
        IH = imread(fileName);
        [r,c,ch] = size(IH);
        
        I3 = imresize(I3, [r,c]);
        IH(:,:,3) = I3;
        
        imwrite(IH, fileName, 'png');
    end
    
end

