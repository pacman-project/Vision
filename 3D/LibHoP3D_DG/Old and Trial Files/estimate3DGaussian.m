% this function is to estimate the parameters of 3D gaussian
% mu is a cnter, Sigma is a covariance matrix

function [mu, Sigma] = estimate3DGaussian(slice)
    
    % slice must be of the odd size
    [xm, ym, zm] = size(slice);
    sumAll = sum(sum(sum(slice)));
    mu = [0, 0, 0];
    halfSize = floor([xm,ym,zm]/2);
    mS = ceil([xm,ym,zm]/2); % middle slice
    
    for i = -halfSize(1):halfSize(1)
        for j = -halfSize(2):halfSize(2)
            for k = -halfSize(3):halfSize(3)
                cur = [i,j,k] * slice(i + mS(1), j + mS(2), k + mS(3));
                mu = mu + cur;
            end
        end
    end
    mu = mu / sumAll;
    
    % now we estimate a covariance matrix
    Sigma = zeros(3);
    
    for i = -halfSize(1):halfSize(1)
        for j = -halfSize(2):halfSize(2)
            for k = -halfSize(3):halfSize(3)
                curX = [i,j,k];
                p = curX - mu;
                curMatr = p'*p;
                curMatr = curMatr * slice(i + mS(1), j + mS(2), k + mS(3));
                Sigma = Sigma + curMatr;
            end
        end
    end
    Sigma = Sigma / sumAll;
    
%     if Sigma(1,1) == 0
%         % sumYZ = (Sigma(2,2) + Sigma(3,3))/100; 
%         Sigma(1,1) = 0.001;
%     end
    


end

