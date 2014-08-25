% this function performs mean value reconstruction of the part of each
% layer

% fieldCenter - position of the central subpart

function [ positions, elements ] = partMeanReconstruction(layerID, partID, fieldCenter, triple8OutDepth, triple7OutDepth, triple6OutDepth, triple5OutDepth, triple4OutDepth, ...
                                                            triple3OutDepth, displ3, displ5, displ7, nClusters)
    
    offsetY4 = [-displ3, 0, displ3];
    offsetX5 = [-displ5, 0, displ5];
    offsetY6 = [-displ5, 0, displ5];
    offsetX7 = [-displ7, 0, displ7];
    offsetY8 = [-displ7, 0, displ7];
    
    indsEls4 = [1,5,6];

    if layerID == 2
        
        positions = fieldCenter;
        elements = partID;
        
    elseif layerID == 3
        
        curEl = triple3OutDepth(partID,:);
        indsEl = [1,2,6,7,8,9];
        indsDepths = [3,4,5,10,11,12];
        curXY = curEl(indsEl);
        depths = curEl(indsDepths);
        
        elements = zeros(1,3);
        for j = 1:3
            elements(j) = compute2elementIndex(curXY(2*j-1), curXY(2*j), nClusters);  % part index in range [0, n2Clusters]
        end
        % define positions
        positionLeft = [fieldCenter(1) - displ3, fieldCenter(2), fieldCenter(3) + depths(3)];
        positionRight =[fieldCenter(1) + displ3, fieldCenter(2), fieldCenter(3) + depths(6)];
        positions = [positionLeft; fieldCenter; positionRight]; % left, middle, right
        
    elseif layerID == 4
        elements = [];
        positions = [];
        
        cur4triple = triple4OutDepth(partID, :);
        els4 = cur4triple(indsEls4);
        depths4Adder = [cur4triple(4), 0 cur4triple(9)];  % top centre bottom
        
        for i = 1:3                            % three triples of the third layet
            curEl = triple3OutDepth(els4(i),:);
            indsEl = [1,2,6,7,8,9];
            indsDepths = [3,4,5,10,11,12];

            curXY = curEl(indsEl);
            depths = curEl(indsDepths);

            for jj = 1:3
                elements = [elements, compute2elementIndex(curXY(2*jj-1), curXY(2*jj), nClusters)];  % part index in range [0, n2Clusters]
            end
            
            positionLeft   = [fieldCenter(1) - displ3, fieldCenter(2) + offsetY4(i), fieldCenter(3) + depths(3) + depths4Adder(i)];
            positionCentre = [fieldCenter(1)         , fieldCenter(2) + offsetY4(i), fieldCenter(3)             + depths4Adder(i)];
            positionRight =  [fieldCenter(1) + displ3, fieldCenter(2) + offsetY4(i), fieldCenter(3) + depths(6) + depths4Adder(i)];

            positions = [positions; positionLeft; positionCentre; positionRight]; % left, middle, right
        end
        
    elseif layerID == 5
        
        elements = [];
        positions = [];
        
        cur5triple = triple5OutDepth(partID, :);
        els5 = cur5triple(indsEls4);
        depths5Adder = [cur5triple(4), 0 cur5triple(9)];  % top centre bottom
        
        
        for k = 1:3 % three triples of the 4th layer
        
            cur4triple = triple4OutDepth(els5(k), :);
            els4 = cur4triple(indsEls4);
            depths4Adder = [cur4triple(4), 0 cur4triple(9)];  % top centre bottom
        
            for i = 1:3                            % three triples of the third layer
                curEl = triple3OutDepth(els4(i),:);
                indsEl = [1,2,6,7,8,9];
                indsDepths = [3,4,5,10,11,12];

                curXY = curEl(indsEl);
                depths = curEl(indsDepths);
                
                for jj = 1:3
                    elements = [elements, compute2elementIndex(curXY(2*jj-1), curXY(2*jj), nClusters)];  % part index in range [0, n2Clusters]
                end

                positionLeft   = [fieldCenter(1) - displ3 + offsetX5(k), fieldCenter(2) + offsetY4(i), fieldCenter(3) + depths(3) + depths4Adder(i) + depths5Adder(k)];
                positionCentre = [fieldCenter(1)          + offsetX5(k), fieldCenter(2) + offsetY4(i), fieldCenter(3)             + depths4Adder(i) + depths5Adder(k)];
                positionRight =  [fieldCenter(1) + displ3 + offsetX5(k), fieldCenter(2) + offsetY4(i), fieldCenter(3) + depths(6) + depths4Adder(i) + depths5Adder(k)];

                positions = [positions; positionLeft; positionCentre; positionRight]; % left, middle, right

            end
            
        end
    
    elseif layerID == 6
        
        elements = [];
        positions = [];
        
        cur6triple = triple6OutDepth(partID, :);
        els6 = cur6triple(indsEls4);
        depths6Adder = [cur6triple(4), 0 cur6triple(9)];  % top centre bottom
        
        
        for kk = 1:3 % three triples of the layer 5
            
            cur5triple = triple5OutDepth(els6(kk), :);
            els5 = cur5triple(indsEls4);
            depths5Adder = [cur5triple(4), 0 cur5triple(9)];  % top centre bottom
        
            for k = 1:3 % three triples of the 4th layer

                cur4triple = triple4OutDepth(els5(k), :);
                els4 = cur4triple(indsEls4);
                depths4Adder = [cur4triple(4), 0 cur4triple(9)];  % top centre bottom

                for i = 1:3                            % three triples of the third layer
                    curEl = triple3OutDepth(els4(i),:);
                    indsEl = [1,2,6,7,8,9];
                    indsDepths = [3,4,5,10,11,12];

                    curXY = curEl(indsEl);
                    depths = curEl(indsDepths);

                    for jj = 1:3
                        elements = [elements, compute2elementIndex(curXY(2*jj-1), curXY(2*jj), nClusters)];  % part index in range [0, n2Clusters]
                    end

                    positionLeft   = [fieldCenter(1) - displ3 + offsetX5(k), fieldCenter(2) + offsetY4(i) + offsetY6(kk), fieldCenter(3) + depths(3) + depths4Adder(i) + depths5Adder(k) + depths6Adder(kk)];
                    positionCentre = [fieldCenter(1)          + offsetX5(k), fieldCenter(2) + offsetY4(i) + offsetY6(kk), fieldCenter(3)             + depths4Adder(i) + depths5Adder(k) + depths6Adder(kk)];
                    positionRight =  [fieldCenter(1) + displ3 + offsetX5(k), fieldCenter(2) + offsetY4(i) + offsetY6(kk), fieldCenter(3) + depths(6) + depths4Adder(i) + depths5Adder(k) + depths6Adder(kk)];

                    positions = [positions; positionLeft; positionCentre; positionRight]; % left, middle, right

                end

            end
        end
        
    elseif layerID == 7
        
        elements = [];
        positions = [];
        
        cur7triple = triple7OutDepth(partID, :);
        els7 = cur7triple(indsEls4);
        depths7Adder = [cur7triple(4), 0 cur7triple(9)];  % top centre bottom
        
        for kkk = 1:3 % three triples of the layer 6
            
            cur6triple = triple6OutDepth(els7(kkk), :);
            els6 = cur6triple(indsEls4);
            depths6Adder = [cur6triple(4), 0 cur6triple(9)];  % top centre bottom
        
            for kk = 1:3 % three triples of the layer 5

                cur5triple = triple5OutDepth(els6(kk), :);
                els5 = cur5triple(indsEls4);
                depths5Adder = [cur5triple(4), 0 cur5triple(9)];  % top centre bottom

                for k = 1:3 % three triples of the 4th layer

                    cur4triple = triple4OutDepth(els5(k), :);
                    els4 = cur4triple(indsEls4);
                    depths4Adder = [cur4triple(4), 0 cur4triple(9)];  % top centre bottom

                    for i = 1:3                            % three triples of the third layer
                        curEl = triple3OutDepth(els4(i),:);
                        indsEl = [1,2,6,7,8,9];
                        indsDepths = [3,4,5,10,11,12];

                        curXY = curEl(indsEl);
                        depths = curEl(indsDepths);

                        for jj = 1:3
                            elements = [elements, compute2elementIndex(curXY(2*jj-1), curXY(2*jj), nClusters)];  % part index in range [0, n2Clusters]
                        end

                        positionLeft   = [fieldCenter(1) - displ3 + offsetX5(k) + offsetX7(kkk), fieldCenter(2) + offsetY4(i) + offsetY6(kk), fieldCenter(3) + depths(3) + depths4Adder(i) + depths5Adder(k) + depths6Adder(kk) + depths7Adder(kk)];
                        positionCentre = [fieldCenter(1)          + offsetX5(k) + offsetX7(kkk), fieldCenter(2) + offsetY4(i) + offsetY6(kk), fieldCenter(3)             + depths4Adder(i) + depths5Adder(k) + depths6Adder(kk) + depths7Adder(kk)];
                        positionRight =  [fieldCenter(1) + displ3 + offsetX5(k) + offsetX7(kkk), fieldCenter(2) + offsetY4(i) + offsetY6(kk), fieldCenter(3) + depths(6) + depths4Adder(i) + depths5Adder(k) + depths6Adder(kk) + depths7Adder(kk)];

                        positions = [positions; positionLeft; positionCentre; positionRight]; % left, middle, right

                    end

                end
            end
        end
        
        elseif layerID == 8
        
        elements = [];
        positions = [];
        
        cur8triple = triple8OutDepth(partID, :);
        els8 = cur8triple(indsEls4);
        depths8Adder = [cur8triple(4), 0 cur8triple(9)];  % top centre bottom
        
        for kkkk = 1:3 % three triples of the layer 7

            cur7triple = triple7OutDepth(els8(kkkk), :);
            els7 = cur7triple(indsEls4);
            depths7Adder = [cur7triple(4), 0 cur7triple(9)];  % top centre bottom
        
            for kkk = 1:3 % three triples of the layer 6

                cur6triple = triple6OutDepth(els7(kkk), :);
                els6 = cur6triple(indsEls4);
                depths6Adder = [cur6triple(4), 0 cur6triple(9)];  % top centre bottom

                for kk = 1:3 % three triples of the layer 5

                    cur5triple = triple5OutDepth(els6(kk), :);
                    els5 = cur5triple(indsEls4);
                    depths5Adder = [cur5triple(4), 0 cur5triple(9)];  % top centre bottom

                    for k = 1:3 % three triples of the 4th layer

                        cur4triple = triple4OutDepth(els5(k), :);
                        els4 = cur4triple(indsEls4);
                        depths4Adder = [cur4triple(4), 0 cur4triple(9)];  % top centre bottom

                        for i = 1:3                            % three triples of the third layer
                            curEl = triple3OutDepth(els4(i),:);
                            indsEl = [1,2,6,7,8,9];
                            indsDepths = [3,4,5,10,11,12];

                            curXY = curEl(indsEl);
                            depths = curEl(indsDepths);

                            for jj = 1:3
                                elements = [elements, compute2elementIndex(curXY(2*jj-1), curXY(2*jj), nClusters)];  % part index in range [0, n2Clusters]
                            end

                            positionLeft   = [fieldCenter(1) - displ3 + offsetX5(k) + offsetX7(kkk), fieldCenter(2) + offsetY4(i) + offsetY6(kk), fieldCenter(3) + depths(3) + depths4Adder(i) + depths5Adder(k) + depths6Adder(kk) + depths7Adder(kkk) + depths8Adder(kkkk)];
                            positionCentre = [fieldCenter(1)          + offsetX5(k) + offsetX7(kkk), fieldCenter(2) + offsetY4(i) + offsetY6(kk), fieldCenter(3)             + depths4Adder(i) + depths5Adder(k) + depths6Adder(kk) + depths7Adder(kkk) + depths8Adder(kkkk)];
                            positionRight =  [fieldCenter(1) + displ3 + offsetX5(k) + offsetX7(kkk), fieldCenter(2) + offsetY4(i) + offsetY6(kk), fieldCenter(3) + depths(6) + depths4Adder(i) + depths5Adder(k) + depths6Adder(kk) + depths7Adder(kkk) + depths8Adder(kkkk)];

                            positions = [positions; positionLeft; positionCentre; positionRight]; % left, middle, right

                        end

                    end
                end
            end  
        end
        
    end
        
    


end

