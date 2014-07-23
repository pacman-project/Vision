function upsampleImages(list_depth, list_mask, lenF, is_downsampling, dowsample_rate, dataSetNumber)

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

            [r, c] = size(I);
            
%             if dataSetNumber == 2
%                 maxDim = max(r,c);
%                 dowsample_rate = 200/maxDim;
%             end


            if ismask
                mask = imread(list_mask{i});
            end

            I = imresize(I, dowsample_rate);
            if ismask
                mask = imresize(mask, dowsample_rate);
            end           
                        
            [r,c] = size(I);

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