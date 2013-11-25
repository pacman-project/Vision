function compute_2D_to_2_5D_correlations()

close all;clc;

image_folder = '/home/c7031079/matlabProjects/3D_lhop/TrimedImages';
images_2D_files = '/home/c7031079/LHOP/Examples/ETH80/testImages/vis_res';
images_2_5_D = '/home/c7031079/matlabProjects/3D_lhop/Elements2Layer';

folder_im = dir([image_folder '/' '*.png']);

%parameters
treshold = 0.01;
borderSize = 120;
delta = 10;
nrLayer_2D = 3;
scales_2D = [0 1 2];
patchSize = 3;
nrParts_2D = [6 20 25 87 342];
nrParts_2_5D_layer3 = 160;

%correlation matrix
correlations_2D = zeros(nrParts_2D(nrLayer_2D),nrParts_2_5D_layer3);

for i=1:size(folder_im)
    
    imName = folder_im(i).name;
    imageName = imName(1:end-4);

    im_2_5D = imread([images_2_5_D '/' imName ]);
    [xpos_2_5D,ypos_2_5D,value_2_5D] = find(im_2_5D');
    im_2_5D_ = im_2_5D';
    
    im_2D_text = [images_2D_files '/test_' imageName '_' num2str(nrLayer_2D) '_' num2str(scales_2D(1)) '.txt'  ]; 
    f = fopen(im_2D_text);
    
    if (f>0)
          
      mat = textscan(f,'%f,%f,%f,%f,%f,%f,%f\n');
      matrix = cell2mat(mat);
      
      for j=1:size(matrix,1)
          
          xpos = matrix(j,1);
          ypos = matrix(j,2);
          type_part2D = matrix(j,4);
          
          %correct coordinates by taking into account borders and scaling
          xpos_ = round((xpos - borderSize)*0.85+delta+delta/2);
          ypos_ = round((ypos - borderSize)*0.85+delta+delta/2);
          
         if((xpos_>0)&&(ypos_>0))
         
          xpos = ypos_; 
          ypos = xpos_;
           
          %rectangle around position in 2D 
          typePart_2_5D = findCorrelation(xpos,ypos,patchSize,im_2_5D_);
          if (typePart_2_5D >0)
              correlations_2D(type_part2D, typePart_2_5D) = correlations_2D(type_part2D, typePart_2_5D) + 1;
          end;
         end;    
      end;
     fclose(f);
     end;
end;

nr_parts_2D = 0;
parts2D = zeros(nrParts_2D((nrLayer_2D)),1);

for i=1:size(correlations_2D,1)
   sum_ = sum(correlations_2D(i,:)); 
   if(sum_>0)
    correlations_2D(i,:) = correlations_2D(i,:)./sum_;
    parts2D(i) = 1; 
    nr_parts_2D = nr_parts_2D + 1;
   end;
end;

%compute log likelihood of correlations greater than treshold
correlations_2D(correlations_2D<treshold)=0;
vectMat = reshape(correlations_2D,1,[]);
norm(log(vectMat(vectMat>0)))

end

function typePart_2_5D = findCorrelation(xpos,ypos,patchSize,im_2_5D_)

typePart_2_5D=0;

if(xpos-patchSize>0) xbegin = xpos-patchSize;
elseif (xpos>0) xbegin = xpos;
else xbegin=0;
end;
         
if(ypos-patchSize>0) ybegin = ypos-patchSize;
elseif (ypos>0) ybegin = ypos;
else ybegin=0;
end;    

if(xpos+patchSize<=size(im_2_5D_,1)) xend = xpos+patchSize;
elseif (xpos<=size(im_2_5D_,1)) xend = xpos;
else xend=0;
end;
         
if(ypos+patchSize<=size(im_2_5D_,2)) yend = ypos+patchSize;
elseif (ypos<=size(im_2_5D_,2)) yend = ypos;
else yend=0;
end;    
         
%build a rectangle around position in 2D 
if(xbegin&&xend&&ybegin&&yend)
 rect = im_2_5D_(xbegin:xend,ybegin:yend);
 [x,y] = find(rect);
 if (size(x)>0)
     typePart_2_5D=im_2_5D_(xbegin+x(1)-1,ybegin+y(1)-1);
 end;
end;

end         