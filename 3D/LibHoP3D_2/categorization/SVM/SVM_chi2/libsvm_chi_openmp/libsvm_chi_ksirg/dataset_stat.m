%data set statistics
clc
dataFolder = 'D:\UWM\doktorat\code\KMLib\KMLibUsageApp\Data\';

tr_path='kytea-msr.train';
tst_path='kytea-msr.test';

tr_path=[dataFolder tr_path];
tst_path=[dataFolder tst_path];


[trYY trXX]=libsvmread(tr_path);
[tstYY tstXX]=libsvmread(tst_path);

%% train stats
DS=trXX;

[~, cols]=find(DS);
maxIdx = max(cols);

sparsity=nnz(DS)/numel(DS);

nnzSum =full( sum(DS' ~=0));

avgEl =mean( nnzSum);
stdEl =std(nnzSum);

[maxNNZ idx1] = max(nnzSum);
[minNNZ idx2] = min( nnzSum);
numRows=size(DS);
fprintf('dataset %s\n stats\n rows=%d dim=%d \n ',tr_path,numRows(1),numRows(2));
fprintf('sparsity=%g \n mean nnz El=%g\n std nnz El=%g',sparsity,avgEl,stdEl);
fprintf('\n min El=%d \n max El=%d\n max nnz idx=%d \n',minNNZ,maxNNZ,maxIdx);


%% test stats
DS=tstXX;

[~, cols]=find(DS);
maxIdx = max(cols);

sparsity=nnz(DS)/numel(DS);

nnzSum =full( sum(DS' ~=0));

avgEl =mean( nnzSum);
stdEl =std(nnzSum);

[maxNNZ idx1] = max(nnzSum);
[minNNZ idx2] = min( nnzSum);

numRows=size(DS);
fprintf('dataset %s\n stats\n rows=%d dim=%d \n ',tst_path,numRows(1),numRows(2));
fprintf('sparsity=%g \n mean nnz El=%g\n std nnz El=%g',sparsity,avgEl,stdEl);
fprintf('\n min El=%d \n max El=%d\n max nnz idx=%d \n',minNNZ,maxNNZ,maxIdx);


