% this function computes areas of the triangles

function [areas] = compute_Areas(V, F)
    
    F1= F(1,:); F2 = F(2,:); F3 = F(3,:);
    V1 = V(:,F1); V2 = V(:,F2); V3 = V(:,F3);
    
    a = sqrt(sum((V1-V2).^2));
    b = sqrt(sum((V1-V3).^2));
    c = sqrt(sum((V2-V3).^2));
    p = (a + b + c) / 2;
    
    areas = sqrt((p-a).*(p-b).*(p-c).*p);

end

% function dist = vectDist(V1, V2)
%     dist = sqrt(sum((V1-V2).^2));
% end
% 
% %     lenF = size(F, 2);
% %     areas = zeros(1, lenF);
% %     
% %     
% %     parfor i = 1:lenF
% %         a = vectDist(V(:,F(1,i)), V(:, F(2,i)));
% %         b = vectDist(V(:,F(1,i)), V(:, F(3,i)));
% %         c = vectDist(V(:,F(2,i)), V(:, F(3,i)));
% %         
% %         p = (a + b + c) / 2;
% %         areas(i) = sqrt((p-a)*(p-b)*(p-c)*p); 
% %     end



