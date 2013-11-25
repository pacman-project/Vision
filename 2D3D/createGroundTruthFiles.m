function createGroundTruthFiles()

close all;clear all;clc;

objectName = 'mug';
folder = ['/home/c7031079/LHOP/Examples/ETH80/' objectName '_1/train'];

images = dir([folder '/' '*.png']);

 for i=1:size(images)
    
        im = imread([folder '/' images(i).name]); 
       
        %save object region to file
        imageName = images(i).name(1:end-4);
        textFile = fopen([folder '/' imageName '_' objectName '.groundtruth'],'w');
        objectRegion = [5,5,size(im,2)-5,size(im,1)-5];
        
        fprintf(textFile,'%f %f %f %f %s\n',objectRegion,objectName);
        
        fclose(textFile);
 end;

end