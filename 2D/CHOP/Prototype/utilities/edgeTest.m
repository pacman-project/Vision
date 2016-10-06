function [ numberOfDiscrepancies ] = edgeTest( graphLevel )
     adjInfoArr = {graphLevel.adjInfo};
     numberOfDiscrepancies = 0;
     for graphLevelItr = 1:numel(graphLevel)
        curAdjacentNodes = adjInfoArr{graphLevelItr};
        if isempty(curAdjacentNodes)
             continue;
        end
        curAdjacentNodes = curAdjacentNodes(:,2);
        validEdgeIdx = ones(numel(curAdjacentNodes),1) > 0;
        for secItr = 1:numel(validEdgeIdx)
             tempArr = adjInfoArr{curAdjacentNodes(secItr)};
             if ~isempty(tempArr)
                  tempArr = tempArr(:,2);
             end
             if ~ismembc(int32(graphLevelItr), tempArr)
                 validEdgeIdx(secItr) = 0;
             end
        end
        if nnz(validEdgeIdx) ~= numel(validEdgeIdx)
             numberOfDiscrepancies = numberOfDiscrepancies + nnz(~validEdgeIdx);
        end
     end
end

