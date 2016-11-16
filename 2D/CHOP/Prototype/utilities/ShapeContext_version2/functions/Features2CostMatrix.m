function C=Features2CostMatrix(F1,F2,Points1,Points2,maxdist)
F1=F1'; F2=F2';
C=zeros([size(F1,1) size(F2,1)]);
for i=1:size(F1,2)
    P=bsxfun(@plus,F1(:,i)+eps,F2(:,i)');
    M=bsxfun(@minus,F1(:,i),F2(:,i)');
    C=C+((M.^2)./P);
end
P1=bsxfun(@minus,Points1(:,1),Points2(:,1)').^2;
P2=bsxfun(@minus,Points1(:,2),Points2(:,2)').^2;
P=P1+P2;
C(P>maxdist^2)=inf;


