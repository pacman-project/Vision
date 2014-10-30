% Read images and labels.
images=loadMNISTImages('train-images-idx3-ubyte')';
labels = loadMNISTLabels('train-labels-idx1-ubyte');

trainingImagePerClass = 500;
testImagePerClass = 25;

for itr = 0:9
    if ~exist([pwd '/vocab/' num2str(itr)], 'dir')
       mkdir([pwd '/vocab/' num2str(itr)]);
    end
    if ~exist([pwd '/test/' num2str(itr)], 'dir')
       mkdir([pwd '/test/' num2str(itr)]);
    end
end

%Training data.
allSampleIdx = [];
allTestSampleIdx = [];
for itr = 0:9
    sampleIdx = find(labels == itr);
    sampleIdx = datasample(sampleIdx, trainingImagePerClass+testImagePerClass, 'Replace', false);
    allSampleIdx = [allSampleIdx; sampleIdx(1:trainingImagePerClass)];
    allTestSampleIdx = [allTestSampleIdx; sampleIdx((trainingImagePerClass+1):end)];
end

for sampleItr = allSampleIdx'
    img = reshape(images(sampleItr,:), [28, 28]);
    imwrite(img, [pwd '/vocab/' num2str(labels(sampleItr)) '/' num2str(sampleItr) '.png']);
end

%images=loadMNISTImages('t10k-images-idx3-ubyte')';
%labels = loadMNISTLabels('t10k-labels-idx1-ubyte');

for sampleItr = allTestSampleIdx'
    img = reshape(images(sampleItr,:), [28, 28]);
    imwrite(img, [pwd '/test/' num2str(labels(sampleItr)) '/' num2str(sampleItr) '.png']);
end
