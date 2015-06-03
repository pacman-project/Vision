% this function computes curvatures and estimates the Darboux frame form depth images

function marksOut = computeCurvatureDI(I, mask, options, zeroThresh)
    
    elementRadius = options.elementRadius;
    numScales = length(elementRadius);
    
    
    marksOut = zeros(size(I));
    
    [sizeY, sizeX, ch] = size(I);
    if ch > 1
        I = I(:,:,1);
    end
    
    vectColors = eye(3);
    

    for x = elementRadius(1) + 1: options.step: sizeX - elementRadius(1)
       for y = elementRadius(1) + 1: options.step: sizeY - elementRadius(1)   % for each point of the image
           
            vecLen = 10;

            if mask(y,x) == 0 || I(y,x) == 0
                continue;
            end
            
            curDepth = I(y,x);
            
            for i = 1:numScales
                % find the largest scale at which the path can be considered as
                % a planar surface

                [indsXOut, indsYOut, depths] = computeNeighbors(I, x, y, elementRadius(i), options.is2D, options.minValue, options.maxValue);
                
                scatter3(indsXOut, indsYOut, depths);
                hold on

                if length(indsXOut) < options.ignoreThresh * (pi*elementRadius(i)^2);
                    continue;
                end
                

                % now we convert these coordinates to the local frame of
                % reference (centered at the point [j,i,curDepth])
                xs = indsXOut - x;
                ys = indsYOut - y;
                zs = depths - I(y,x);

                if options.methodId == 1  % paraboloid fitting
                    
                    [absErr] = PlanarTestParaboloid(xs, ys, zs);
%                     [H, K, QF1, QF2, S, V, D] = paraboloidFitting(xs, ys, zs);



                elseif options.methodId == 2 % method of Gabriel Taubin (modified)
                    
                    disp('Sorry Gabriel Taubins method is not implemented');
%                      [H, K, V, D] = computeDarbouxFrame(Ix(y,x), Iy(y,x), xs, ys, zs);
                end
            end
 
            if  marksOut(y,x) ~= 0
                
%                 % plot a circle in 3D
%                 plotCircle3D([x,y,curDepth], V(:,1)', marksOut(y,x));
%                 a = 2;
%                 
%                 vecLen = vecLen * (2 * marksOut(y,x) + 1);

                vecLen = 40;
                % plot the eigenvectors in 3D
                for ii = 1:3
                    curVect = V(:, ii);
                    curColor = vectColors(ii, :);
                    XX = [x, x + vecLen * curVect(1)];
                    YY = [y, y + vecLen * curVect(2)];
                    ZZ = [curDepth, curDepth + vecLen * curVect(3)];

                    plot3(XX, YY, ZZ, 'Color', curColor);
                    hold on
                    

                end
            end


            
        a = 2;

        end 
    end
    
end

