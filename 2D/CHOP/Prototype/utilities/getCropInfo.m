function [ cropBounds ] = getCropInfo( exportArr, categoryArr, trainingFileNames )
    imageIds = exportArr(:,5);
    numberOfImages = max(imageIds);
    for imageItr = 1:numberOfImages
       fileName = trainingFileNames{imageItr};
       positions = exportArr(imageIds==imageItr,2:3);
       mu = mean(positions);
    end
end

