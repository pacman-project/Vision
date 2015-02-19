% make all images 3-dimensional

function checkImages(list_depth, lenF)

    disp('making input images 3-dimensional...');
    
    parfor i = 1:lenF
        I = imread(list_depth{i}); 
        I(:,:,2) = I(:,:,1);
        I(:,:,3) = I(:,:,1)+1;   % no mask in this channel. This helps for selection of parts with empty cells at image borders
        imwrite(I, list_depth{i}, 'png');
    end
end