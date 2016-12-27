function [selectedParts, reconstructiveParts, discriminativeParts] = selectJointParts(newCliques, newMetaData, leafNodes, level1Coords, categoryArrIdx, options)
     % Program variables set here.
     discriminativeParts = [];
     coverageThr = 0.95;
     numberOfSubs = double(max(newMetaData(:,1)));
     numberOfImages = double(numel(categoryArrIdx));
     numberOfSelectedParts = options.numberOfSelectedSubs;
     
     %% Create the structures for future processing.
     uniqueMetaData = unique(newMetaData(:,1:2), 'rows');
     [~, idx] = unique(uniqueMetaData(:,1), 'R2012a');
     idx = cat(1, idx, size(uniqueMetaData,1)+1);
     d = cell(numberOfSubs,1);
     parfor itr = 1:numberOfSubs
          dummyArr = zeros(numberOfImages, 1) > 0;
          dummyArr(uniqueMetaData(idx(itr):(idx(itr+1)-1),2)) = 1;
          d{itr} = sparse(dummyArr);
     end
     d = cat(2, d{:});
     
     %% Select parts based on coverage.
     % Calculate remaining nodes per image.
     remainingNodes = newCliques(newCliques>0);
     remainingNodes = fastsortedunique(sort(remainingNodes(:)));
     remainingLeafNodes = fastsortedunique(sort(cat(2, leafNodes{remainingNodes})));
     imageNodeCounts = hist(level1Coords(remainingLeafNodes,1), 1:numberOfImages);
     minCoverages = imageNodeCounts * coverageThr;
     
     % First, save support of parts.
     coveredNodes = cell(numberOfSubs,1);
     valArr = zeros(numberOfSubs, 1);
     [~, idx] = unique(newMetaData(:,1), 'R2012a');
     idx = cat(1, idx, size(newMetaData,1)+1);
     parfor partItr = 1:numberOfSubs
          dummyArr = zeros(1, size(level1Coords,1))>0;
          children = newCliques(idx(partItr):(idx(partItr+1)-1),:);
          children = children(children > 0);
          children = children(:);
          dummyArr(cat(2, leafNodes{children})) = 1;
          nodes = find(dummyArr);
          valArr(partItr) = numel(nodes);
          coveredNodes{partItr} = nodes;
     end
       
     % Start selecting.
     validNodes = ones(max(remainingLeafNodes),1) > 0;
     [~, reconstructiveParts] = max(valArr);
     valArr(reconstructiveParts) = 0;
     validNodes(coveredNodes{reconstructiveParts}) = 0;
     if options.reconstructivePartSelection
          while numel(reconstructiveParts) < numberOfSelectedParts
               % Update valid nodes if some images are covered.
               newImageNodeCounts = hist(level1Coords(~validNodes,1), 1:numberOfImages);
               coveredImages = newImageNodeCounts >= minCoverages & minCoverages > 0;
               validNodes(ismembc(level1Coords(:,1), int32(find(coveredImages)))) = 0;

               % Remove ineligible subs.
               invalidSubs = ~full(sum(d(~coveredImages, :),1)) > 0;
               valArr(invalidSubs) = 0;

               % Go through the list of parts.
               [sortedValArr, sortIdx] = sort(valArr, 'descend');
               eligibleSubs = find(sortedValArr > 0);
               maxVal = 0;
               maxPart = 0;
               for partItr = 1:numel(eligibleSubs)
                    % If we have no hope of getting a better part, close this
                    % loop.
                    if maxVal > sortedValArr(partItr)
                         break;
                    end

                    % Calculate value and save if necessary.
                    tempArr = coveredNodes{sortIdx(partItr)};
                    tempArr = tempArr(validNodes(tempArr));
                    coveredNodes{sortIdx(partItr)} = tempArr;
                    newVal = numel(tempArr);
                    valArr(sortIdx(partItr)) = newVal;

                    % Save the part, it's our champion!
                    if newVal > maxVal
                         maxValidNodes = validNodes;
                         maxValidNodes(tempArr) = 0;
                         maxPart = sortIdx(partItr);
                         maxVal = newVal;
                    end
               end

               % If we haven't found a new part, stop processing.
               if ~maxVal
                    break;
               else
                    reconstructiveParts = cat(1, reconstructiveParts, maxPart);
                    validNodes = maxValidNodes;
               end
          end
          reconstructiveParts = sort(reconstructiveParts);
     else
          reconstructiveParts = [];
     end
     
     if options.discriminativePartSelection
          %% Start Discriminative selection.
          d = double(d);
          nd = size(d,2);

          t1=cputime;
          t = zeros(nd, 1);
          parfor i=1:nd, 
             t(i) = mutualinfo(full(d(:,i)), categoryArrIdx);
          end; 
          fprintf('calculate the marginal dmi costs %5.1fs.\n', cputime-t1);

          [~, idxs] = sort(-t);

          fea = zeros(min(numberOfSelectedParts, nd),1);

          if numberOfSelectedParts < nd
              fea(1) = idxs(1);

              KMAX = min(20000,nd); %500 %20000

              if KMAX <= numberOfSelectedParts
                  fea = idxs((1:numberOfSelectedParts));
                  return;
              end

              idxleft = idxs(2:KMAX);

              k=1;
              fprintf('k=1 cost_time=(N/A) cur_fea=%d #left_cand=%d\n', ...
                    fea(k), length(idxleft));

              mi_array = zeros(max(idxleft), numberOfSelectedParts - 1);
              valArr = zeros(numberOfSelectedParts-1, 1);

              for k=2:numberOfSelectedParts,
                 t1=cputime;
                 ncand = length(idxleft);
                 curlastfea = nnz(fea);
                 t_mi = t(idxleft);
                 temp_array = zeros(ncand, 1);

                 % We calculate mutual information for promising subs.
                 comparedFeature = full(d(:,fea(curlastfea)));
                 parfor i=1:ncand,
                      temp_array(i) = getmultimi(comparedFeature, full(d(:,idxleft(i))));
                 end
                 mi_array(idxleft, curlastfea) = temp_array;
                 c_mi = sum(mi_array(idxleft, 1:(k-1)), 2) / (k-1);

                 % Calculate value and pick the sub with maximum value.
                 [valArr(k-1), fea(k)] = max(t_mi(1:ncand) ./ (c_mi(1:ncand) + 0.01));
                 tmpidx = fea(k); 
                 fea(k) = idxleft(tmpidx); 
                 idxleft(tmpidx) = [];

                 if rem(k, 100) == 0
                 fprintf('k=%d cost_time=%5.4f cur_fea=%d #left_cand=%d\n', ...
                    k, cputime-t1, fea(k), length(idxleft));
                 end
              end

          else
              fea = (1:nd)';
          end
          discriminativeParts = fea;
          discriminativeParts = sort(discriminativeParts);
     else
          discriminativeParts = [];
     end
     %% Combine discriminative parts + reconstructive parts.
     selectedParts = sort(unique(cat(1, reconstructiveParts, discriminativeParts)));
end

%===================================== 
function c = getmultimi(da, dt) 
     for i=1:size(da,2), 
        c(i) = mutualinfo(da(:,i), dt);
     end; 
end