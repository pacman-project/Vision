


dataFolder = 'D:\UWM\doktorat\code\KMLib\KMLibUsageApp\Data\';

tr_path=[dataFolder 'a9a.small.t'];
tst_path=[dataFolder 'a9a.t'];

[trYY trXX]=libsvmread(tr_path);
[tstYY tstXX]=libsvmread(tst_path);

maxCore = 8;
maxIter= floor(log2(maxCore))+1;

cores=ones(1,maxIter);
timeTest=zeros(1, maxIter);
timeTrain = zeros(1, maxIter);

for k=1:maxIter
    
    nrThreads = 2^(k-1);
    % new 'z' - option, determine the number of threads 
    optStr = sprintf('-c 4 -t 2 -g 0.5 -q -z %d',nrThreads);
    tic;
    model = svmtrain(trYY, trXX,optStr);
    modelTime=toc;
    
    tic
    [pred, acc, dec_vals] = svmpredict(tstYY, tstXX, model);
    predTime = toc;
    
    cores(k)=nrThreads;
    timeTrain(k)=modelTime;
    timeTest(k) = predTime;
    ss=sprintf('libsvm RBF acc=%0.5g modeltime=%g predtime=%g \n\n',acc(1),modelTime, predTime);
    disp(ss);
end

figure(1);
hold on
plot(cores,timeTrain,'ro:');
plot(cores, timeTest,'bs-');
hold off;



