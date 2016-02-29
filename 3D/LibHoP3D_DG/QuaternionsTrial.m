% this is another trial with quaternions

function QuaternionsTrial()
    
    numIter = 60;
    numTrain = 40;
    frameCananical = eye(3);  % right-handed system
    
    
    % for visualization of the frame
    n = 3;
    vectColors = eye(3);
    vecLen = 0.3;
    centre = [0,0,0];
    lim = 1.0;
    
%     plotFrame(frameCananical, centre, 3* vecLen, vectColors, n, 1);
%     xlim([-lim lim])
%     ylim([-lim lim])
%     zlim([-lim lim])
%     hold on
    
    qAll = zeros(numIter, 4);
    VAll = zeros(numIter, 3);
    UAll = zeros(numIter, 3);
    NAll = zeros(numIter, 3);
    
    
    % generate a set of right-handed frames with pre-defined Z axis
    
    orientCanonical = [1,1,1];
    vectorStandard = [0,0,1];
    frameStandard = eye(3);
    
    for i = 1:numIter
        
        alpha = 2 * pi * rand(1);
        
        if i <= numTrain
            
            V = orientCanonical + 0.5 * rand(1,3);
            V = V/norm(V);
            
            % define a quaternion
            q = computeQuaternion(vectorStandard, V);
            q = qnorm(q);
            qAll(i, :) = q;
            continue;
            

% %             V = [cos(alpha), sin(alpha), 0];
% %             U = [sin(alpha), cos(alpha), 0];
% %             N = [0,0,1];
% % 
% %             temp = dot(cross(V,U), N);% check the right-handiness
% % 
% %             if ~(0.99 < temp && temp < 1.01)
% %                 V = -V;
% %             end
% 
%             V = [1,1,1] + 0.2 * rand(1,3);
%             V = V/norm(V);
% %             Utemp = zeros(1,3);
% %             r =find(abs(V) == min(abs(V)));
% %             Utemp(r) = 1;
%             Utemp = [0,0,1];
%             N = cross(V, Utemp);
%             N = N/norm(N);
%             U = cross(N, V);
%             U = U/norm(U);

        else
            
            N = orientCanonical + 0.2 * rand(1,3);
            N = N/norm(N);
            Vtemp = [1,0,0];
            U = cross(N, Vtemp);
            U = U/norm(U);
            V = cross(U, N);
            V = V/norm(V);
            
            a= 2;
            
%             V = orientCanonical - 2 * rand(1,3);
%             V = V/norm(V);
%             Utemp = zeros(1,3);
%             r =find(abs(V) == min(abs(V)));
%             Utemp(r) = 1;
%             N = cross(V, Utemp);
%             N = N/norm(N);
%             U = cross(N, V);
%             U = U/norm(U);

        end
        
        
%         VAll(i, :) = V;
%         UAll(i, :) = U;
%         NAll(i, :) = N;
%         
%         plotFrame([V', U', N'], centre, vecLen, vectColors, n, 1);
%         hold on
        
        % compute directional cosine matrix and quaternion
        dcMatr = [V;U;N];
        q = dcm2q(dcMatr);
        q = qnorm(q);
        
%         addRand = 0 * (1 - 2*rand(1,4))/100;  % range [-0.03, 0.03]
%         q = q + addRand;
%         q = qnorm(q);
        
        qAll(i, :) = q;
    end

    GMModel = fitgmdist(qAll(1:numTrain, :), 1, 'RegularizationValue', 10^-13);
    mu = GMModel.mu;
    Sigma = GMModel.Sigma;
    invSigma = inv(Sigma);
    
%     scoreGauss(i,jj) = sqrt((points(jj, 1:3) - curMu) * invSigma * (points(jj, 1:3) - curMu)');
    score = zeros(1, numIter);
    for i = 1:numIter
        score(i) = sqrt((qAll(i, :) - mu) * invSigma * (qAll(i, :) - mu)');
    end
    
    sum(score(1:numTrain))/length(score(1:numTrain))
    sum(score(1+numTrain:end))/length(score(1+numTrain:end))
    a = 2;
    
    
    
%     for i = 1:numIter
%         % generate a random (right-handed) frame of reference
%         V = 1 - 2*rand(1,3);
%         V = V/norm(V);
% 
%         Utemp = zeros(1,3);
%         r =find(abs(V) == min(abs(V)));
%         Utemp(r) = 1;
% 
%         N = cross(V, Utemp);
%         N = N/norm(N);
%         U = cross(N, V);
%         U = U/norm(U);
%         
%         VAll(i, :) = V;
%         UAll(i, :) = U;
%         NAll(i, :) = N;
% % 
%         plotFrame([V', U', N'], centre, vecLen, vectColors, n, 1);
%         hold on
% 
%         % compute directional cosine matrix and quaternion
%         dcMatr = [V;U;N];
%         q = dcm2q(dcMatr);
%         q = qnorm(q);
%         qAll(i, :) = q;
%     end
    
%     % check smoothness of quaternions
%     for i = 1:numIter
%         
%         figure
%         q = qAll(i, :);
%         
%         addRand = 3 * (1 - 2*rand(1,4))/100;  % range [-0.03, 0.03]
%         q = q + addRand;
%         q = qnorm(q);
% %         
% %         % original frame of reference
%         V = VAll(i, :);
%         U = UAll(i, :);
%         N = NAll(i, :);
%         
%         Vq = qvrot(q, frameCananical(:, 1))';
%         Uq = qvrot(q, frameCananical(:, 2))';
%         Nq = qvrot(q, frameCananical(:, 3))'
%         
%         plotFrame([V', U', N'], centre, vecLen, vectColors, n, 1);
%         hold on
%         plotFrame([Vq', Uq', Nq'], centre, vecLen/2, vectColors, n, 2);
%         
%         a = 2;
%         
%         % transform the original frame of reference using quaternion
%         
%         
%     end
    
    
    
    % clustering of the quaternions
%     nCl = 20;
%     cl = kmeans(qAll,nCl);
%     figure
%     
%     for i = 1:nCl
%         ids = find(cl == i);
% %         figure
%         numvects = length(ids);
%         
%         vectColors = rand(3,3);
%         for j = 1:numvects
%             V = VAll(ids(j), :);
%             U = UAll(ids(j), :);
%             N = NAll(ids(j), :);
%             plotFrame([V',U', N'], centre, vecLen, vectColors, n, 1);
%             
%             xlim([-lim lim])
%             ylim([-lim lim])
%             zlim([-lim lim])
%             hold on
%         end
%         a = 2;
%     end

end

function plotFrame(V, centre, vecLen, vectColors, n, lineWidth)

    % plot the eigenvectors in 3D
    for ii = 1:n
        curVect = V(:, ii);
        curColor = vectColors(ii, :);
        XX = [centre(1), centre(1) + vecLen * curVect(1)];
        YY = [centre(2), centre(2) + vecLen * curVect(2)];
        ZZ = [centre(3), centre(3) + vecLen * curVect(3)];

        plot3(XX, YY, ZZ, 'Color', curColor, 'LineWidth', lineWidth);
        hold on

    end
end