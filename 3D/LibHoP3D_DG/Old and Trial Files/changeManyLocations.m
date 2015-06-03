

function [] = changeManyLocations()
    
    % initialize imagesTrial
    numImages = 1800;
    imSize = 600;
    dx = 7;
    dy = 3;
    
    tic

    parfor i = 1:numImages
        I = uint8(zeros(imSize, imSize));
        filename = ['ImagesTrial/', num2str(i), '.png'];
        
        imwrite(I, filename, 'png');
    end
    
    len = 100;

%     imIds = gpuArray.randi(numImages, numImages*len, 1);
%     x = gpuArray.randi(500, numImages*len, 1) + 50;
%     y = gpuArray.randi(500, numImages*len, 1) + 50;
    
    imIds = randi(numImages, numImages*len, 1);
    x = randi(500, numImages*len, 1) + 50;
    y =randi(500, numImages*len, 1) + 50;
    
    toc
    
    pointFun = @givePiece;
    
    tic
    
    
    parfor i = 1:numImages
        
        inds = find(imIds == i);
        
        if ~isempty(inds)

            filename = ['ImagesTrial/', num2str(i), '.png'];
            I = imread(filename);

            for j = 1:length(inds)
                I(y(inds(j))-dy:y(inds(j))+dy, x(inds(j))-dx:x(inds(j))+dx) = 255;  
            end

            I = uint8(I);
            imwrite(I, filename, 'png');
        end
        
    end

    
%     parfor i = 1:numImages
%         
%         [xxx, yyy] = feval(pointFun, i);
%         
%         if ~isempty(xxx)
% 
%             filename = ['ImagesTrial/', num2str(i), '.png'];
%             I = imread(filename);
% 
%             for j = 1:length(xxx)
%                 I(yyy(j)-dy:yyy(j)+dy, xxx(j)-dx:xxx(j)+dx) = 255;  
%             end
% 
%             I = uint8(I);
%             imwrite(I, filename, 'png');
%         end
%         
%     end
%     
%     
%     function [xx, yy] = givePiece(i)
%         indSS = imIds == i;
%         xx = gather(x(indSS));
%         yy = gather(y(indSS));  
%     end
    
    
    toc
    
    


end

