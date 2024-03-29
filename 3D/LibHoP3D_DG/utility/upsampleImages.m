function upsampleImages(list_depth, list_mask, lenF, is_downsampling, dowsample_rate)

    disp('upsampling images...');
    
    if ~isempty(list_mask)
        ismask = true;
    else
        ismask = false;
        list_mask = zeros(1,lenF);
    end
    
    if is_downsampling
        
        parfor i = 1:lenF
            I = imread(list_depth{i});
            if ismask
                mask = imread(list_mask{i});
            end
            
            I = double(I) * dowsample_rate; % to preserve geometry
            I = imresize(I, dowsample_rate);
            if ismask
                mask = imresize(mask, dowsample_rate);
            end                        
            
            I = uint16(I);
            imwrite(I, list_depth{i}, 'png');
            if ismask
                imwrite(mask, list_mask{i}, 'png');
            end

            if mod(i,100) == 0
                i
            end
        end

    end
end