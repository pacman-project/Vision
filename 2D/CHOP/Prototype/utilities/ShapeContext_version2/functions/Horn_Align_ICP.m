function [P2,P2L] = Horn_Align_ICP(P1,P2,P2L)
if(nargin<3), P2L=P2; end

check=true; Cost=inf;
while(check);
    Costold=Cost;
    ns1 = createns(P1,'nsmethod','kdtree');
    ns2 = createns(P2,'nsmethod','kdtree');
    M1 = knnsearch(ns1,P2,'k',1);
    M2 = knnsearch(ns2,P1,'k',1);
    Vertices=[P2;P2(M2,:)];
    VerticesB=[P1(M1,:);P1];
    Cost=mean((sum((Vertices-VerticesB).^2,2)));
    disp(['Horn Point Distance: ' num2str(Cost)]); drawnow;
    [P2,P2L]=Horn_Align_ICP2D3D(Vertices,VerticesB,P2,P2L);
    check=Costold>Cost*1.02;  
end


function [P2N,P2LN] =Horn_Align_ICP2D3D(Vertices,VerticesB,P2,P2L)
is3D=size(Vertices,2)==3;
% Center data to remove translation 
offsetv = -mean(Vertices,1);
Vertices(:,1) = Vertices(:,1) + offsetv(1);
Vertices(:,2) = Vertices(:,2) + offsetv(2);
if(is3D)
    Vertices(:,3) = Vertices(:,3) + offsetv(3);
end

offsetvB = -mean(VerticesB,1);
VerticesB(:,1) = VerticesB(:,1) + offsetvB(1);
VerticesB(:,2) = VerticesB(:,2) + offsetvB(2);
if(is3D)
    VerticesB(:,3) = VerticesB(:,3) + offsetvB(3);
end

% Warp-Matrix between point sets
S=Vertices'*VerticesB;
if(is3D)
    M=[(S(1,1)+S(2,2)+S(3,3))  (S(2,3)-S(3,2))      (S(3,1)-S(1,3))      (S(1,2)-S(2,1));...
       (S(2,3)-S(3,2))      (S(1,1)-S(2,2)-S(3,3))  (S(1,2)+S(2,1))      (S(3,1)+S(1,3));...
       (S(3,1)-S(1,3))      (S(1,2)+S(2,1))     (-S(1,1)+S(2,2)-S(3,3))  (S(2,3)+S(3,2));...
       (S(1,2)-S(2,1))      (S(3,1)+S(1,3))      (S(2,3)+S(3,2))      (-S(1,1)-S(2,2)+S(3,3))];
else
    M=[(S(1,1)+S(2,2)) (S(2,1)-S(1,2));...
        (S(2,1)-S(1,2)) -(S(1,1)+S(2,2))];
end

%Compute eigenvalues
[E,~] = eig(M);

% Calculate the quaternion
if(is3D)
    quat = E(:,4);
else
    quat = E(:,2);    
end
[~,ind]=max(abs(quat)); 
quat=quat.*sign(quat(ind)); 
quat=quat./norm(quat);

%Compute the rotation matrix
if(is3D)
    q0=quat(1); qx=quat(2);  qy=quat(3); qz=quat(4); qxyz =quat(2:4);
    R=qxyz*qxyz.' + [q0 -qz qy; qz q0 -qx; -qy qx  q0]^2;
else
    R=[(quat(1)^2-quat(2)^2) (quat(1)*quat(2)*2);-(quat(1)*quat(2)*2) (quat(1)^2-quat(2)^2)]; 
end
VerticesR=(R*Vertices')';
summ = @(M) sum(M(:));
offsets=summ( VerticesB.*VerticesR)/summ(VerticesR.^2);
R=R*offsets;
T = -offsetvB' + R*offsetv';

% Make transformation matrix
M2=[R T]';

% Warp the vertices
if(is3D)
    M2(4,4)=1;
    P2(:,4)=1; P2N=P2*M2; P2N=P2N(:,1:3);
    P2L(:,4)=1; P2LN=P2L*M2; P2LN=P2LN(:,1:3);
else
    M2(3,3)=1;
    P2(:,3)=1; P2N=P2*M2; P2N=P2N(:,1:2);
    P2L(:,3)=1; P2LN=P2L*M2; P2LN=P2LN(:,1:2);
end
