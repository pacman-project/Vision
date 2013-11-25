function display_2D_to_2_5D_correlations()

close all;clc;

image_folder = '/home/c7031079/matlabProjects/3D_lhop/TrimedImages';
images_2D_save = '/home/c7031079/matlabProjects/3D_lhop/2D_images_layer4';
images_2_5_D = '/home/c7031079/matlabProjects/3D_lhop/Elements2Layer';
fused_images = '/home/c7031079/matlabProjects/3D_lhop/fusedImages';

folder_im = dir([image_folder '/' '*.png']);
patchSize = 7;

for i=1:size(folder_im)
    
    imageName = folder_im(i).name;
    im = imread([image_folder '/' imageName]);
    
    im_2_5D = imread([images_2_5_D '/' imageName ]);
    [xpos_2_5D,ypos_2_5D,value_2_5D] = find(im_2_5D');
    
    figure;imshow(im);
    for j=1:size(value_2_5D)
        hold on;
        circleplot(xpos_2_5D(j),ypos_2_5D(j),2,1); 
    end;
    %hold off;
    
    if exist([images_2D_save '/' imageName ], 'file') == 2

     %figure;imshow(im);
     im_2D = imread([images_2D_save '/' imageName ]);
     [xpos_2D,ypos_2D,value_2D] = find(im_2D');
 
     for j=1:size(value_2D)
        hold on;
        circleplot(xpos_2D(j),ypos_2D(j),2,2); 
     end;
     %hold off;
    
     %figure; imshow(im);
     
     %%display correlations from 2D -> 2.5D
     im_2_5D_ = im_2_5D';
     for j=1:size(value_2D)
        
        %find correlation around coordinates in 2D
        if (findCorrelation(xpos_2D(j),ypos_2D(j),patchSize,im_2_5D_))
           hold on
           circleplot(xpos_2D(j),ypos_2D(j),2,3);
        end;
        
     end;
     
     %%display correlations from 2_5D -> 2D
     im_2D_ = im_2D';
     for j=1:size(value_2_5D)
       
       %find correlation around coordinates in 2_5D   
       if (findCorrelation(xpos_2_5D(j),ypos_2_5D(j),patchSize,im_2D_))
           hold on
           circleplot(xpos_2_5D(j),ypos_2_5D(j),2,4);
       end;
        
     end;
     
    end;
end; 
    
end

function result = findCorrelation(xpos,ypos,patchSize,im_2_5D_)

result=false;

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
     result = true;
 end;
end;

end         
         
function circleplot(xc, yc, r,color)

t = 0 : .1 : 2*pi;
x = r * cos(t) + xc;
y = r * sin(t) + yc;

if(color==1)
  
    plot(x, y,'y');
  
elseif (color==2)
    
  plot(x, y,'m');

elseif (color==3)

    plot(x, y,'b');

elseif (color==4)

    plot(x, y,'g');  

end;

axis square; grid 

end