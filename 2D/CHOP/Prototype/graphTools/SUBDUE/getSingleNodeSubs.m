%> Name: getSingleNodeSubs
%>
%> Description: getSingleNodeSubs is used to obtain single-node subs from
%> labels of all instances in allLabels. The result is a number of
%> subs representing compositions each having their instances.
%> 
%> @param allLabels Labels for every graph node.
%> @param allSigns Signs for every graph node.
%> @param categoryArrIdx Categories for every graph node.
%> @param validationIdx Binary array showing if each node is part 
%> of the validation data or not.
%>
%> @retval singleNodeSubs Substructure list of single-node subs.
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 24.02.2014
%> Ver 1.1 on 01.09.2014 Removal of global parameters.
function singleNodeSubs = getSingleNodeSubs(allLabels, allSigns, nodeDistanceMatrix, categoryArrIdx, validationIdx, threshold)
    numberOfSubs = max(allLabels);
    singleNodeSubs(numberOfSubs) = Substructure();
    validSubs = ones(numberOfSubs,1)>0;
    
    %% For each center node label type, we create a substructure.
    for subItr = 1:numberOfSubs
        distances = nodeDistanceMatrix(subItr, allLabels)';
        subCenterIdx = distances < threshold;
        instances = int32(find(subCenterIdx));
        numberOfInstances = numel(instances);
        
        %Assign center id.
        singleNodeSubs(subItr).centerId = subItr;
        
        % Give maximum score so that it is at the top of the queue.
        singleNodeSubs(subItr).mdlScore = numberOfInstances;
        
        if numberOfInstances>0
            % Fill in instance information. 
            categoryIdx = categoryArrIdx(subCenterIdx);
            if size(categoryIdx,1) == 1
                categoryIdx = categoryIdx';
            end
            instanceValidationIdx = validationIdx(subCenterIdx);
            instanceSigns = allSigns(subCenterIdx,1);
            instanceDistances = distances(subCenterIdx);
            instanceMappings = ones(numberOfInstances, 1, 'uint8');
            
            % We check if this sub has exact-matching instances in the
            % training set. If not, it is not considered for further expansion.
            if nnz(instanceDistances == 0) == 0
                 singleNodeSubs(subItr).mdlScore = 0;
                 continue;
            end
            
            % Assign fields of the sub.
            singleNodeSubs(subItr).instanceCenterIdx = instances;
            singleNodeSubs(subItr).instanceChildren = instances;
            singleNodeSubs(subItr).instanceMappings = instanceMappings;
            singleNodeSubs(subItr).instanceCategories = categoryIdx;
            singleNodeSubs(subItr).instanceMatchCosts = instanceDistances;
            singleNodeSubs(subItr).instanceExactMatchFlags = allLabels(instances) == subItr;
            singleNodeSubs(subItr).instanceSigns = instanceSigns;
            singleNodeSubs(subItr).instanceValidationIdx = instanceValidationIdx;
        else
            validSubs(subItr) = 0;
        end
    end
    
    % Eliminate those with no instances.
    singleNodeSubs = singleNodeSubs(validSubs);
end