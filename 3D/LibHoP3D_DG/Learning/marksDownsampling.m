% this function performs downsampling of MARKS images before learning of each layer
% downSamplingFactor can be 2 or 3

function [] = marksDownsampling(list_El, lenF, downSamplingFactor, elPathIn, elPathOut)

    disp('Downsampling of images with parts ... ');
    
% lenF = 1;
% downSamplingFactor = 3;
% list_El{1} = 'C:\Projects\Vladislav\Input data\Mirela_dataset\Mirela_dataset_layer4\objects_depth\bowl_ceramic_white\bowl_ceramic_white_16.png';

    lenDPW = length(elPathIn);

    parfor i = 1:lenF   
        
        marks = imread(list_El{i});
        [r,c] = size(marks);
        
        
        if downSamplingFactor == 2
            
            outMarks = zeros(floor(r/2), floor(c/2));
            
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
                end
            end
            
        elseif downSamplingFactor == 3  % this is supposed to be applied after layer 4
            
            outMarks = zeros(floor(r/3), floor(c/3));
            
            for ii = 2:3:r
                for jj = 2:3:c
                    
%                     if jj+1<=c && ii+1 <= r
%                         wind = marks(ii-1:ii+1, jj-1:jj+1);
%                     else
%                         wind = marks(ii-1:r, jj-1:c);
%                     end

                    wind = marks(ii-1:min(r, ii+1), jj-1:min(c, jj+1));
                    
                    a = wind(:);
                    a(a == 0) = [];
                    if isempty(a)
                        mm = 0;
                    else
                        mm = mode(a);
                    end
                    outMarks((ii+1)/3, (jj+1)/3) = mm;
                end
            end
            
            
        end
        
        if mod(i,10) == 0
            i
        end

%         imtool(marks, [min(min(marks)), max(max(marks))]);
%         imtool(outMarks, [min(min(outMarks)), max(max(outMarks))]);
        
        curStr = list_El{i};
        fileName = curStr(lenDPW+1:end);
        fileOut = [elPathOut, fileName]; 
         
        % check if this folder exists
        ll = strfind(fileOut, '/');
        lll = strfind(fileOut, '\');
        ll = [ll, lll];
        ll = max(ll); % last position
        folderName = fileOut(1:ll-1);
  
        b = exist(folderName,'dir');
        if b == 0
            mkdir(folderName);
        end
        
        outMarks = uint16(outMarks);
        imwrite(outMarks, fileOut, 'png');
     
    end

end

