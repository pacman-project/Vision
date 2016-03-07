% This is to compute distance from each point to each cluster centre (Mahalanobis distance + Vladislav's distance)

function D = Mixed_Eucl_Quat_cluster_dist(points, parts, alpha)

    if nargin == 2
        alpha = 0.75;
    end
    
    numParts = size(parts, 1);
    numPointsT = size(points, 1);
    scoreGauss = zeros(numParts, numPointsT);
    scoreFischer = zeros(numParts, numPointsT);
    
    scoreGauss1 = zeros(numParts, numPointsT);
    scoreFischer1 = zeros(numParts, numPointsT);
    nCl = size(parts, 1);

    for i = 1:nCl

        curMu = parts(i, 1:3); % centre of each clusters
        curMuOrient = parts(i, 4:7);
        curSigma = reshape(parts(i, 8:16), [3,3]);
        curSigmaOrient = reshape(parts(i, 17:32), [4,4]);
        invSigma = inv(curSigma);
        invSigmaOrient = inv(curSigmaOrient);
        
% %         if numPointsT > 10^5
% %             parfor jj = 1:numPointsT
% %                 scoreGauss(i,jj) = sqrt((points(jj, 1:3) - curMu) * invSigma * (points(jj, 1:3) - curMu)');
% % %                 scoreFischer(i,jj) = 1 - (vonMisesMultiplier * exp(kappa*points(jj, 4:6)*mu')/kappa) / vonMisesDenominator; % Vladislav's distance
% %                 scoreFischer(i,jj) = sqrt((points(jj, 4:7) - curMuOrient) * invSigmaOrient * (points(jj, 4:7) - curMuOrient)');
% %             end
% %         else
% %             for jj = 1:numPointsT
% %                 scoreGauss(i,jj) = sqrt((points(jj, 1:3) - curMu) * invSigma * (points(jj, 1:3) - curMu)');
% % %                 scoreFischer(i,jj) = 1 - (vonMisesMultiplier * exp(kappa*points(jj, 4:6)*mu')/kappa) / vonMisesDenominator; % Vladislav's distance
% %                 scoreFischer(i,jj) = sqrt((points(jj, 4:7) - curMuOrient) * invSigmaOrient * (points(jj, 4:7) - curMuOrient)');
% %             end
% %         end
        
        tXYZ = bsxfun(@minus, points(:, 1:3), curMu);
        tQ = bsxfun(@minus, points(:, 4:7), curMuOrient);
        scoreGauss(i, :) =  sqrt(sum(tXYZ * invSigma .* tXYZ, 2))';
        scoreFischer(i, :) = sqrt(sum(tQ * invSigmaOrient .* tQ, 2))';
        
    end
    D = (scoreGauss + alpha * scoreFischer)';
    a = 2;

%     scoreAll = scoreGauss + 3 *scoreFischer;
%     [m,cEst] = min(scoreAll(1:3, :), [], 1);
%     sc = [scoreAll; scoreGauss; scoreFischer; c'; cEst; m];  % for testing
end









% % %     for i = 1:nCl
% % %         curMu = parts(i, 1:3);
% % %         kappa = parts(i, 16);
% % %         mu = parts(i, 4:6);
% % %         curSigma = reshape(parts(i, 7:15), [3,3]);
% % %         invSigma = inv(curSigma);
% % %         vonMisesMultiplier = kappa/(2 * pi * (exp(kappa) - exp(-kappa)));
% % %         vonMisesDenominator = vonMisesMultiplier * exp(kappa*mu*mu')/kappa;
% % %         
% % %         if numPointsT > 10^5
% % %             parfor jj = 1:numPointsT
% % %                 scoreGauss(i,jj) = sqrt((points(jj, 1:3) - curMu) * invSigma * (points(jj, 1:3) - curMu)');
% % %                 scoreFischer(i,jj) = 1 - (vonMisesMultiplier * exp(kappa*points(jj, 4:6)*mu')/kappa) / vonMisesDenominator; % Vladislav's distance
% % %             end
% % %         else
% % %             for jj = 1:numPointsT
% % %                 scoreGauss(i,jj) = sqrt((points(jj, 1:3) - curMu) * invSigma * (points(jj, 1:3) - curMu)');
% % %                 scoreFischer(i,jj) = 1 - (vonMisesMultiplier * exp(kappa*points(jj, 4:6)*mu')/kappa) / vonMisesDenominator; % Vladislav's distance
% % %             end
% % %         end
% % %     end
% % %     D = (scoreGauss + alpha * scoreFischer)';

