%> Name: evaluateSubs
%>
%> Description: The evaluation function of SUBDUE. Based on the
%> evaluation metric, the value of each substructure in subs list is
%> calculated, and saved within subs. The MDL calculation takes place here.
%> Key points:
%>  1) DL estimation not rigorous. It is considered that 
%>        each node needs a label (int) and a pointer to its edges (int)
%>        each edge needs a node label (int) for its destination node, an
%>        edge label(int) and a binary isDirected label (bit)
%>     in the final graph description. Node and edge weights in DL
%>     calculation is stored in options.
%> 
%> @param subs Sub list which will be evaluated.
%> @param evalMetric Evaluation metric.
%> @param allEdges List of all edges in the graph.
%> @param allEdgeNodePairs List of all edge node pairs in the graph.
%> @param allSigns List of all signs of the nodes in the graph.
%> @param graphSize Size of the graph.
%> @param overlap If true, overlapping instances are considered in
%> evaluation of a sub.
%> @param mdlNodeWeight Node weight in DL calculations (MDL).
%> @param mdlEdgeWeight Edge weight in DL calculations (MDL).
%> @param isMDLExact If true, exact MDL calculation. Approximate otherwise.
%>
%> @retval subs Evaluated sub list.
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 24.02.2014
%> Ver 1.1 on 01.09.2014 Removal of global parameters.
function [subs] = evaluateSubs(subs, evalMetric, allEdges, allEdgeNodePairs, allSigns, graphSize, overlap, mdlNodeWeight, mdlEdgeWeight, isMDLExact, isSupervised, isMDLNormalized)
    numberOfSubs = numel(subs);
    parfor subItr = 1:numberOfSubs
        % Find the weight of this node, by taking the max of the category distribution. 
        if isSupervised
            categoryArr = double(subs(subItr).instanceCategories);
            weight = nnz(categoryArr == mode(categoryArr)) / numel(categoryArr);
        else
            weight = 1;
        end
        
        % We compress the object graph using the children, and the
        % edges they are involved. 
        [subScore, sub, numberOfNonoverlappingInstances] = getSubScore(subs(subItr), allEdges, allEdgeNodePairs, evalMetric, ...
           allSigns, mdlNodeWeight, mdlEdgeWeight, ....
            overlap, isMDLExact);
        subScore = subScore * weight;
        subs(subItr) = sub;
        
        %% Assign the score of the sub, as well as its normalized mdl score if applicable.
        subs(subItr).mdlScore = double(subScore);
        if strcmp(evalMetric, 'mdl') || strcmp(evalMetric, 'likelihood')
            subs(subItr).normMdlScore = 1 - (double(subScore) / graphSize);
            
            % If the normalized MDL score is being asked for, we divide
            % subScore by the number of valid instances.
            if isMDLNormalized
                 subs(subItr).mdlScore = subs(subItr).mdlScore / numberOfNonoverlappingInstances;
            end
        end
    end
end