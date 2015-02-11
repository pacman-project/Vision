% this function performs downsampling of MARKS images before learning of each layer
% downSamplingFactor can be 2 or 3

function [ ] = marksDownsampling(list_depths, list_El, list_mask, lenF, dataSetNumber, downSamplingFactor, isErrosion, discSize,...
                                is_mask_extended, maxExtThresh1, maxExtThresh2, isMaskResize)

    disp('Marks downsampling');
    
    for i = 1:lenF   
        
%         I = imread(list_depths{i});
%         I = I(:,:,1);
%       
%         if dataSetNumber == 2
%             
%             mask = imread(list_mask{i});
%             [~, ~, ~, mask, r, c, is_successfull] = preliminaryProcessing(I, mask, isErrosion, discSize, false, false, false, ...
%                 false, [], 0, 0, 2, 0, 0, is_mask_extended, maxExtThresh1, maxExtThresh2, [], [], [], []);
%             
%             
%         else
%             % to define mask
%             
%             [~, ~, ~, mask, r, c, is_successfull] = preliminaryProcessing(I, [], isErrosion, discSize, false, false, false, ...
%                     false, [], 0, 0, 2, 0, 0, false, 0, 0, ...
%                     [],[],[],[]);
%         end
        
        marks = imread(list_El{i});
        
        I = marks + 100;
        I(I == 100) = 0;
        imtool(I, [1, 200]);
        
        [r,c] = size(marks);
%         imtool(marks, [0, 50]);
%         imtool(mask, [0, 1]);
        
        if downSamplingFactor == 2
            
            outMarks = zeros(floor(r/2), floor(c/2));
            if isMaskResize
                outMask = false(floor(r/2), floor(c/2));
            end
            
            for ii = 2:2:r
                for jj = 2:2:c
                    
                    wind = marks(ii-1:ii, jj-1:jj);
                    a = wind(:);
                    a(a == 0) = [];
                    if isempty(a)
                        mm = 0;
                    else
                        mm = mode(a);
                    end
                    outMarks(ii/2, jj/2) = mm;
                    
                    if isMaskResize
                        % resize a mask here
                        wind = mask(ii-1:ii, jj-1:jj);
                        a = wind(:);
                        mm = length(a(a == 1));
                        if mm > 2
                            outMask(ii/2, jj/2) = true;
                        end
                    end
                end
            end
            
        elseif downSamplingFactor == 3  % this is supposed to be applied after layer 4
            
            outMarks = zeros(floor(r/3), floor(c/3));
            outMask = false(floor(r/3), floor(c/3));
            
            for ii = 2:3:r
                for jj = 2:3:c
                    
                    if jj+1<=c && ii+1 <= r
                        wind = marks(ii-1:ii+1, jj-1:jj+1);
                    else
                        wind = marks(ii-1:r, jj-1:c);
                    end
                    
                    a = wind(:);
                    a(a == 0) = [];
                    if isempty(a)
                        mm = 0;
                    else
                        mm = mode(a);
                    end
                    outMarks((ii+1)/3, (jj+1)/3) = mm;
                    
                    % resize a mask here
                    if isMaskResize
                        wind = mask(ii-1:ii, jj-1:jj);
                        a = wind(:);
                        mm = length(a(a == 1));
                        if mm > 2
                            outMask((ii+1)/3, (jj+1)/3) = true;
                        end
                    end
                end
            end
            
            
        end
        
        if mod(i,10) == 0
            i
        end
        
%         imtool(outMarks, [0, 50]);
%         imtool(outMask, [0, 1]);
        
        % save the downsampled images in place
        
        I = outMarks + 100;
        I(I == 100) = 0;
        imtool(I, [1, 200]);
        outMarks = uint16(outMarks);
        
        imwrite(outMarks, list_El{i}, 'png');
     
    end

    % TODO something smart with mask

end

