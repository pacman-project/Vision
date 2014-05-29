function [F] = EDGEPS(f)
    % (c) 2010 Nikos Papamarkos
    % Democritus University of Thrace, Greece.
    % papamark@ee.duth.gr
    % This program performs Edge preserving smoothing in color images.
    % The program is based on the algorithm described in the following paper:  
    % N. Nikolaou and N. Papamarkos  "Color reduction for complex document 
    % images", International Journal of Imaging Systems and Technology, 
    % Vol. 19 (1) pp. 14-26, 2009. 
    figure (1), imshow(f), title('Initial Image');
    f=double(f);
    % Usually the number of iteration is 1-10
    iter=2;
    for p=0:iter
        r=f(:,:,1); 
        g=f(:,:,2);
        b=f(:,:,3);
        [m n]=size(r);
    for i=1:m
        for j=1:n
            ip=i+1;im=i-1;jm=j-1;jp=j+1;
            if(im<1)
                im=i;
            elseif (ip>m)
                ip=i;
            end
            if(jm<1)
                jm=j;
            elseif (jp>n)
                jp=j;
            end
            d(1)=(abs(r(i,j)-r(i,jm))+abs(g(i,j)-g(i,jm))+abs(b(i,j)-b(i,jm)));
            d(2)=(abs(r(i,j)-r(i,jp))+abs(g(i,j)-g(i,jp))+abs(b(i,j)-b(i,jp)));
            d(3)=(abs(r(i,j)-r(ip,j))+abs(g(i,j)-g(ip,j))+abs(b(i,j)-b(ip,j)));
            d(4)=(abs(r(i,j)-r(im,j))+abs(g(i,j)-g(im,j))+abs(b(i,j)-b(im,j)));
            d(5)=(abs(r(i,j)-r(im,jm))+abs(g(i,j)-g(im,jm))+abs(b(i,j)-b(im,jm)));
            d(6)=(abs(r(i,j)-r(ip,jm))+abs(g(i,j)-g(ip,jm))+abs(b(i,j)-b(ip,jm)));      
            d(7)=(abs(r(i,j)-r(ip,jp))+abs(g(i,j)-g(ip,jp))+abs(b(i,j)-b(ip,jp)));      
            d(8)=(abs(r(i,j)-r(im,jp))+abs(g(i,j)-g(im,jp))+abs(b(i,j)-b(im,jp)));  
            d=d/(3*255);
            sum2=0;
            for k=1:8
              c(k)=(1-d(k))^p;sum2=sum2+c(k);
            end
            rs(i,j)=floor((1/sum2)*(c(1)*r(i,jm)+c(2)*r(i,jp)+c(3)*r(ip,j)+c(4)*r(im,j)+c(5)*r(im,jm)+c(6)*r(ip,jm)+c(7)*r(ip,jp)+c(8)*r(im,jp)));
            gs(i,j)=((1/sum2)*(c(1)*g(i,jm)+c(2)*g(i,jp)+c(3)*g(ip,j)+c(4)*g(im,j)+c(5)*g(im,jm)+c(6)*g(ip,jm)+c(7)*g(ip,jp)+c(8)*g(im,jp)));
            bs(i,j)=((1/sum2)*(c(1)*b(i,jm)+c(2)*b(i,jp)+c(3)*b(ip,j)+c(4)*b(im,j)+c(5)*b(im,jm)+c(6)*b(ip,jm)+c(7)*b(ip,jp)+c(8)*b(im,jp)));
        end
    end
        F=cat(3,rs,gs,bs);
    end
end