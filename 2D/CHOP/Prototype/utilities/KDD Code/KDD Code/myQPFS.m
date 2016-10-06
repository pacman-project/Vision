%Nguyen X. Vinh, Jeffrey Chan, Simone Romano and James Bailey, "Effective Global Approaches for Mutual Information based Feature Selection". 
%To appear in Proceeedings of the 20th ACM SIGKDD Conference on Knowledge Discovery and Data Mining (KDD'14), August 24-27, New York City, 2014. 
% (C) 2014 Nguyen Xuan Vinh   
% Email: vinh.nguyen@unimelb.edu.au, vinh.nguyenx@gmail.com
% Performing the following mutual informattion based feature selection
% approaches:
% - Maximum relevance (maxRel)
% - Minimum redundancy maximum relevance (MRMR)
% - Minimum redundancy (minRed)
% - Quadratic programming feature selection (QPFS)
% - Mutual information quotient (MIQ)
function [QPFSFS,maxRelFS,MRMRFS,QPFS_MImat,minRedFS,MIQFS,timeMaxRel,timeMRMR,timeQPFS,timeMIQ]=myQPFS(a,C,maxFeature,QPFS_MImat)
[n dim]=size(a);
if nargin<3 maxFeature=dim;end;

H=zeros(dim,dim);
f=zeros(dim,1);
maxE=log(max(max(a)));

fprintf('QPFS started, calculating the MI matrix...')

if nargin<4||isempty(QPFS_MImat) 
    tic;
    H=computeMImatrix_4([a C]);    
    matrixTime=toc;
    QPFS_MImat=H; 
else H=QPFS_MImat;end;

f=H(1:end-1,end);
H=H(1:end-1,1:end-1);

%-doing maxRel and MRMR
tic;[MRMRFS]=mrmr(H,f,maxFeature);timeMRMR=toc()+ matrixTime*(maxFeature/dim ); %MRMR
tic;[minRedFS]=minRed(H,f,maxFeature);timeMaxRed=toc()+ matrixTime*(maxFeature/dim ); %max redundancy
tic;[MIQFS]=MIQ(H,f,maxFeature);timeMIQ=toc()+ matrixTime*(maxFeature/dim ); %Mutual information quotient

%--maxRel
tic;
y=zeros(dim,2);
for i=1:dim y(i,1)=-f(i);y(i,2)=i;end;
y=sortrows(y);
y=[y(:,2) -y(:,1)];
maxRelFS=y(:,1)';
timeMaxRel=toc();

%-doing QPFS
%H=(H+lambda*eye(dim));
tic;
evalue=eig(H);
if evalue(1)<-10^-3
    fprintf('QPFS: non-positive Q\n');
    for i=1:dim H(i,i)=H(i,i)+3;end;
end

mq=sum(sum(H))/dim^2;
mf=sum(f)/dim;
alpha=mq/(mq+mf)
f=-alpha*f;
H=(1-alpha)*H;

A=-eye(dim);  % x_i >=0 contraints
b=zeros(dim,1);
Aeq=ones(1,dim);
beq=1;

options = optimset('Algorithm','interior-point-convex');
fprintf('Solving the QPFS formulation...\n');
x = quadprog(H,f,A,b,Aeq,beq,[],[],[],options);

f=0.5*(1-alpha)*x'*H*x -alpha*f'*x
y=zeros(dim,2);
y(:,1)=-x;
for i=1:dim y(i,2)=i;end;
y=sortrows(y);
y=[y(:,2) -y(:,1)];
QPFSFS=y(:,1);
timeQPFS=toc()+matrixTime;

end


function [MRMRFS]=mrmr(H,f,maxFeature)
[N dim]=size(H);
if nargin<3 maxFeature=dim;end;

fprintf('MRMR started...\n');

%select first feature as the one with max MI with C
max_MI=0;firstFeature=1;
for i=1:dim
    CMI=f(i);
    if CMI>max_MI
        max_MI=CMI;
        firstFeature=i;
    end
end
best_fs=zeros(1,maxFeature);
best_fs(1)=firstFeature;
%fprintf('Adding first feature %d\n',firstFeature);

selected=zeros(1,dim);
selected(best_fs(1))=1;

%forward selection
%fprintf('Forward selection...\n');
for j=2:maxFeature
    max_inc=-inf;
    bestFeature=0;
    for i=1:dim
       if selected(i) continue;end;
       rel=f(i);       
       red=sum(H(i,best_fs(1:j-1)))/(j-1);
       
       inc=rel-red;
       if inc>max_inc
           max_inc=inc;
           bestFeature=i;
       end
    end
    
    best_fs(j)= bestFeature;
    selected(bestFeature)=1;
    %fprintf('Adding %d\n',bestFeature);    
end
MRMRFS=best_fs;
end


function [MRMRFS]=minRed(H,f,maxFeature)  %min redundancy feature selection
[N dim]=size(H);
if nargin<3 maxFeature=dim;end;

fprintf('Min redundancy started...\n');

%select first feature as the one with max MI with C
max_MI=0;firstFeature=1;
for i=1:dim
    CMI=f(i);
    if CMI>max_MI
        max_MI=CMI;
        firstFeature=i;
    end
end
best_fs=zeros(1,maxFeature);
best_fs(1)=firstFeature;
%fprintf('Adding first feature %d\n',firstFeature);

selected=zeros(1,dim);
selected(best_fs(1))=1;

%forward selection
%fprintf('Forward selection...\n');
for j=2:maxFeature
    max_inc=-inf;
    bestFeature=0;
    for i=1:dim
       if selected(i) continue;end;
       rel=0;  %don't care about relevancy       
       red=sum(H(i,best_fs(1:j-1)))/(j-1);
       
       inc=rel-red; %Just min redundancy
       if inc>max_inc
           max_inc=inc;
           bestFeature=i;
       end
    end
    
    best_fs(j)= bestFeature;
    selected(bestFeature)=1;
    %fprintf('Adding %d\n',bestFeature);    
end
MRMRFS=best_fs;
end

function [MRMRFS]=MIQ(H,f,maxFeature)  %mutual information quotient
[N dim]=size(H);
if nargin<3 maxFeature=dim;end;

fprintf('MIQ started...\n');

%select first feature as the one with max MI with C
max_MI=0;firstFeature=1;
for i=1:dim
    CMI=f(i);
    if CMI>max_MI
        max_MI=CMI;
        firstFeature=i;
    end
end
best_fs=zeros(1,maxFeature);
best_fs(1)=firstFeature;
%fprintf('Adding first feature %d\n',firstFeature);

selected=zeros(1,dim);
selected(best_fs(1))=1;

%forward selection
%fprintf('Forward selection...\n');
for j=2:maxFeature
    max_inc=-inf;
    bestFeature=0;
    for i=1:dim
       if selected(i) continue;end;
       rel=f(i);       
       red=sum(H(i,best_fs(1:j-1)))/(j-1);
       
       inc=rel/red;  %MIQ criterion
       if inc>max_inc
           max_inc=inc;
           bestFeature=i;
       end
    end
    %fprintf('Rel= %f red=%f rel/red= %f\n',rel,red,inc);
    best_fs(j)= bestFeature;
    selected(bestFeature)=1;
    fprintf('Adding %d\n',bestFeature);    
end
MRMRFS=best_fs;
end