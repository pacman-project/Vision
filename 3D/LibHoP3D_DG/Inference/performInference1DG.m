% this function computes curvatures and estimates the Darboux frame form depth images

function [marksOut, Normals] = performInference1DG(I, mask, options)
    
    elementRadius = options.elementRadius;
    numScales = length(elementRadius);

    [sizeY, sizeX, ch] = size(I);
    if ch > 1
        I = I(:,:,1);
    end
    
    errors = zeros(sizeY, sizeX);
    marksOut = zeros(sizeY, sizeX);
    Normals = zeros(sizeY, sizeX, 3);

    for x = elementRadius(1) + 1: options.step: sizeX - elementRadius(1)
       for y = elementRadius(1) + 1: options.step: sizeY - elementRadius(1)   % for each point of the image
% 
%             
%             if ~(y == 278 && x== 300)
%                 continue;
%             end
            
            if mask(y,x) == 0 || I(y,x) < options.minValue
                continue;
            end
            
            absErr = 0;
            i = 1;
            is_ok = false;
            
               
            [indsXOut, indsYOut, depths] = computeNeighbors(I, x, y, elementRadius(i), options.is2D, options.minValue, options.maxValue);
            if length(indsXOut) < options.ignoreThresh{1};
                continue;
            end

            % now we convert these coordinates to the local frame of
            % reference (centered at the point [j,i,curDepth])
            xs = indsXOut - x;
            ys = indsYOut - y;
            zs = depths - I(y,x);

            if options.methodId == 1  % paraboloid fitting  
%                         
%               [absErr, Norm] = PlanarTestParaboloid(xs, ys, zs);
                [absErr, Norm, is_ok] = PlanarTestLeastSquares(xs, ys, zs);
                is_ok = true;
                errors(y,x) = absErr;
                
            elseif options.methodId == 2 % method of Gabriel Taubin (modified)
                disp('performInference1DG: Sorry Gabriel Taubins method is not implemented');
            end 

            if (absErr < options.planarPatch) && (i==1)
                marksOut(y,x) = i;
                Normals(y,x,1) = Norm(1); 
                Normals(y,x,2) = Norm(2); 
                Normals(y,x,3) = Norm(3);
            end

       end 
    end
    
%     imtool(Normals);
%     imtool(errors, [0, 0.2]);
%     
%     imtool(marksOut,[0,1]);
%     a=2;
  

end

    
    
    
    
    %                     
%             elseif is_multyScale   % multiple scale representation
%                 
%                 while (absErr < options.planarPatch) && (i <= numScales)
%                     % find the largest scale at which the path can be considered as
%                     % a planar surface
% 
%                     [indsXOut, indsYOut, depths] = computeNeighbors(I, x, y, elementRadius(i), options.is2D, options.minValue, options.maxValue);
% 
%                     if length(indsXOut) < options.ignoreThresh{1} * (pi*elementRadius(i)^2);
%                         i = i + 1;
%                         continue;
%                     end
% 
%                     % now we convert these coordinates to the local frame of
%                     % reference (centered at the point [j,i,curDepth])
%                     xs = indsXOut - x;
%                     ys = indsYOut - y;
%                     zs = depths - I(y,x);
% 
%                     if options.methodId == 1  % paraboloid fitting                  
%                         [absErr, Norm, is_ok] = PlanarTestParaboloid(xs, ys, zs);
%                     elseif options.methodId == 2 % method of Gabriel Taubin (modified)
%                         disp('Sorry Gabriel Taubins method is not implemented');
%                     end
%                     i = i+1;
%                 end
% 
%                 i = i-1;
%                 
%                 if is_ok
%                     % check which exit condition holds 
%                     if (i == numScales) && (absErr < options.planarPatch)
%                         marksOut(y,x) = numScales;
%                         Normals(y,x,1) = Norm(1); Normals(y,x,2) = Norm(2); Normals(y,x,3) = Norm(3);
%                     elseif (absErr > options.planarPatch) && (i>1)
%                         marksOut(y,x) = i;
%                         Normals(y,x,1) = Norm(1); Normals(y,x,2) = Norm(2); Normals(y,x,3) = Norm(3);
%                     elseif (absErr > options.planarPatch) && (i==1)
%                         % no detection at this point
%                     else
%                         disp('Something went wrong in performInference1DG');
%                     end
%                 end
% 
%             end
            


