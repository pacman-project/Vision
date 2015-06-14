% make sure all images are one-dimensional

function checkImages2(list_depth, lenF)

    disp('making input images 1-dimensional...');
    
    parfor i = 1:lenF
        I = imread(list_depth{i}); 
        I = I(:,:,1);
        imwrite(I, list_depth{i}, 'png');
    end
end