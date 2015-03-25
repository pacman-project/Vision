% elementType : 1 - point
%               2 - disk 
%               3 - square

function [learningElType, learningElRadius, inferenceElType, inferenceElRadius] = loadLearningInferenceStructElement(dataSetNumber)

    
    if dataSetNumber == 1 || dataSetNumber == 3

        % learning
        learningElType{3} = 2;
        learningElType{4} = 2;
        learningElType{5} = 2;
        learningElType{6} = 2;
        learningElType{7} = 2;
        learningElType{8} = 2;
        
        learningElRadius{3} = 2;
        learningElRadius{4} = 2;
        learningElRadius{5} = 2;
        learningElRadius{6} = 2;
        learningElRadius{7} = 2;
        learningElRadius{8} = 2;
        
        % inference
        inferenceElType{3} = 2;
        inferenceElType{4} = 2;
        inferenceElType{5} = 2;
        inferenceElType{6} = 2;
        inferenceElType{7} = 2;
        inferenceElType{8} = 2;
        
        inferenceElRadius{3} = 2;
        inferenceElRadius{4} = 2;
        inferenceElRadius{5} = 2;
        inferenceElRadius{6} = 2;
        inferenceElRadius{7} = 2;
        inferenceElRadius{8} = 2;
        
    elseif dataSetNumber == 2
        for i = 3:8
            
            % learning
            learningElType{3} = 2;
            learningElType{4} = 2;
            learningElType{5} = 2;
            learningElType{6} = 2;
            learningElType{7} = 2;
            learningElType{8} = 2;

            learningElRadius{3} = 2;
            learningElRadius{4} = 2;
            learningElRadius{5} = 2;
            learningElRadius{6} = 2;
            learningElRadius{7} = 2;
            learningElRadius{8} = 2;

            % inference
            inferenceElType{3} = 2;
            inferenceElType{4} = 2;
            inferenceElType{5} = 2;
            inferenceElType{6} = 2;
            inferenceElType{7} = 2;
            inferenceElType{8} = 2;

            inferenceElRadius{3} = 2;
            inferenceElRadius{4} = 2;
            inferenceElRadius{5} = 2;
            inferenceElRadius{6} = 2;
            inferenceElRadius{7} = 2;
            inferenceElRadius{8} = 2;
            
        end
        
    end

end
