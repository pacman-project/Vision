% this is the script for different methods of estimating principal
% curvatures

%% image of the sphere to establish relation of depth and x-y units
I = imread('G:\VladislavSTD\Vladislav_STD\depthAll\99_4_2_1_3_6.png');
I = I(:,:,1);
I = double(I);

% determine radius of the sphere
[r,c]= find(I>0);
rMin = min(r);
rMax = max(r);
cMin = min(c);
cMax = max(c);

diam = rMax - rMin;
rad = diam/2;
I= I(rMin:rMax, cMin:cMax);

centerX = round(diam/2);
centerY = round(diam/2);

dDepth = I(centerY, centerX) - I(centerY, centerX + round(diam/4));
scaleD = rad*(2-sqrt(2))/dDepth;
zeroThresh = 10^(-4);

% I = I * scaleD;

% a = 2;
% imtool(I, [min(min(I)), max(max(I))]);

%% Estiamte principal curvatures and darboux frame of the cylinder

I = imread('G:\VladislavSTD\Vladislav_STD\depthAll\99_5_2_1_1_6.png');
I = I(:,:,1);
I = double(I);
I = I * scaleD;

[row,col]= find(I>0);
rMin = min(row);
rMax = max(row);
cMin = min(col);
cMax = max(col);
I= I(rMin:rMax, cMin:cMax); % trim the image

[row,col] = size(I);
methodID = 1;

IH = zeros(row,col);
IK = zeros(row,col);

if methodID == 1  % derivation of all parameters from the analytical surface fitting
    
    % the surface is of type f = z(x,y):
    % z = a*x^2 + b*y^2 + c x y + d x + e y + f
    
    % receptive field size
    radRC = 10;  % in fact the receptive feid is a box of size (2,2,2)
    minValue = 1000;
    
    for i = radRC + 1 : row - radRC
        for j = radRC + 1 : col - radRC
            
            if I(i, j) ~= 0  % check the central element
                
                bBox = I(i-radRC:i+radRC, j-radRC:j+radRC);
                [rr,cc] = find(bBox > minValue);
                % local frame of reference has an origin in the point of the
                % bounding box
                if length(rr) < 100
                    continue;
                end   

                xs = cc - radRC - 1;
                ys = rr - radRC - 1;
                xs2 = xs.^2;
                ys2 = ys.^2;
                xys = xs.*ys;
                ind = sub2ind(size(bBox), rr,cc);
                zs = bBox(ind);

                X = [xs2, ys2, xys, xs, ys, ones(length(rr),1)];
                A = X\zs;  % (a, b, c, d, e, f)

                % estimate the normal in the centre of the receptive field (x = 0, y = 0)
                tangentX = [1,0,A(4)];
                tangentY = [0,1,A(5)];
                Norm = cross(tangentX, tangentY);
                
                Norm = Norm/norm(Norm);
                
                % estimate the first and second quadratic forms
                % we have to compute cos(N,z)
                z = [0,0,1];
                cosNZ = dot(Norm, z);
                QF1 = [1+A(4)^2, A(4)*A(5); A(4)*A(5), 1+A(5)^2];
                QF2 = [2*A(1),2*A(3); 2*A(3), 2*A(2)] * cosNZ;
                
                E = QF1(1,1); F = QF1(2,1); G = QF1(2,2);
                L = QF2(1,1); M = QF2(2,1); N = QF2(2,2);
                
                % gaussian and mean curvatures
                K = (L*M - N^2)/(E*G - F^2);
                H = (E*N + G*L -2*F*M) / (2*(E*G - F^2));
                
                if abs(K) < zeroThresh
                    K = 0;
                end
                if abs(H) < zeroThresh
                    H = 0;
                end
                IH(i,j) = H;
                IK(i,j) = K;
                
                % now we have to compute principal curvatures
                % solve the following quadratic equation a*k^2+b*k + c = 0 where k is unknown
                a = E*G - F^2; 
                b = -(E*N - 2*F*M + G*L); 
                c = L*N - M^2;
                D = b^2 - 2*a*c;  % discriminant
                
                
                if D > zeroThresh % there are two roots
                   k(1) = ( -b + sqrt(D) ) / (2*a);
                   k(2) = ( -b - sqrt(D) ) / (2*a);
                end
                
                % we find the shape operator
                
                S = inv(QF1)*QF2;
                
                [V,D] = eig(S);
                
                if i == 75 && j == 50
                
                    x1 = [j - 10*V(1,1), j, j + 10*V(1,1)];
                    y1 = [j - 10*V(2,1), i, j + 10*V(2,1)];

                    x2 = [j - 10*V(1,2), j, j + 10*V(1,2)];
                    y2 = [j - 10*V(2,2), i, j + 10*V(2,2)];
                    imshow(I);
                    hold on
                    plot(x1,y1, 'blue');
                    hold on
                    plot(x2,y2, 'red'); 
                    
                    a = 2;
                end

            end
        end
    end
    
%     imtool(IH, [min(min(IH)), max(max(IH))]);
%     imtool(IK, [min(min(IK)), max(max(IK))]);
    a = 2;
    
elseif methodID == 1  % here we compute the Dabourx frama directly
    
    
    
end


































