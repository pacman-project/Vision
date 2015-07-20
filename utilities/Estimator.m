%> Name: Estimator
%>
%> Description: The class which calculates likelihood of data points given
%> the reconstructed (backprojected) features. Additionally, it calculates
%> the reconstruction error, i.e. the discrepancy between data points and
%> their corresponding reconstructions. The algorithm assumes that both 
%> dataCoords and reconstructedCoords are of the same size, and they have
%> 1-to-1 (ordered) correspondance. 
%> 
%> Author: rusen
%>
%> Updates
%> Ver 1.0 on 03.07.2015
classdef Estimator
    properties (Access = private)
        numberOfPoints = -1;
        
        % PDF parameters.
        sigma = 0.598;
        epsilon = 0.05;
        concentration = 0.945; % So that Pr(0, pi, concentration) = epsilon.
        
        % Data structures.
        dataCoords@uint16; % Nx2 array of uint16. Each row holds coordinates (x,y) of a single data point.
        dataLabels@uint16; % Nx1 array of uint16. Point labels.
        dataDepths@uint16; % Nx1 array of uint16. Used for 2.5D. Point depths.
        reconstructedCoords@uint16; % Nx2 array of uint16. Each row holds coordinates (x,y) of a reconstructed feature.
        reconstructedLabels@uint16; % Nx1 array of uint16. Reconstructed feature labels.
        reconstructedDepths@uint16; % Nx1 array of uint16. Used for 2.5D. Reconstructed feature depths.
        radius@double;      % Constant. The radius in which the algorithm looks for correspondances.
        numberOfFeatures@uint8; % Constant. Number of low level features.
    end
    
    methods
        %% Core methods
        % The function that estimates data likelihood by comparing data features 
        % and reconstructed features. DataType is a string that can be
        % added to specify different cases (2D, 2.5D, 3D). 
        function likelihood = CalculateDataLikelihood(this)
            likelihoodArr = zeros(this.numberOfPoints, 1);
            
            % Convert labels to angles for von mises distribution.
            anglePerFeature = (2*pi) / double(this.GetNumberOfFeatures());
            dataAngles = double( this.dataLabels) * anglePerFeature;
            reconstructedAngles = double(this.reconstructedLabels) * anglePerFeature;
            
            for itr = 1:this.numberOfPoints
                 % Calculate spatial probability.
                 if ismember(0, this.dataCoords(itr,:)) || ismember(0, this.reconstructedCoords(itr,:))
                      % Invalid coords, which means assigned to an
                      % non-existing node.
                      spatialProb = this.epsilon;
                 else
                      spatialProb = mvnpdf(double(this.dataCoords(itr,:))/this.radius, double(this.reconstructedCoords(itr,:))/this.radius, [this.sigma, this.sigma]);
                 end
                 
                 % Calculate node replacement (appearance) probability.
                 
                 likelihoodArr(itr) = spatialProb * circ_vmpdf(dataAngles(itr,1), reconstructedAngles(itr,1), this.concentration);
            end
            
           % Calculate the final likelihood.
           likelihood = prod(likelihoodArr);
        end
        
        
        % The function that estimates log likelihood by comparing data features 
        % and reconstructed features. DataType is a string that can be
        % added to specify different cases (2D, 2.5D, 3D). 
        function logLikelihood = CalculateLogLikelihood(this)
           logLikelihoodArr = zeros(this.numberOfPoints, 1);
            
            % Convert labels to angles for von mises distribution.
            anglePerFeature = (2*pi) / double(this.GetNumberOfFeatures());
            dataAngles = double( this.dataLabels) * anglePerFeature;
            reconstructedAngles = double(this.reconstructedLabels) * anglePerFeature;
            downsampleRatio = this.radius / sqrt(2);
            for itr = 1:this.numberOfPoints
                 % Calculate spatial probability.
                 if ismember(0, this.dataCoords(itr,:)) || ismember(0, this.reconstructedCoords(itr,:))
                      % Invalid coords, which means assigned to an
                      % non-existing node.
                      spatialProb = log2(this.epsilon);
                      angularProb = spatialProb;
                 else
                      spatialProb = log2(mvnpdf(double(this.dataCoords(itr,:))/downsampleRatio, double(this.reconstructedCoords(itr,:))/downsampleRatio, [this.sigma, this.sigma]));
                      angularProb = log2(circ_vmpdf(dataAngles(itr,1), reconstructedAngles(itr,1), this.concentration));
                 end
                 
                 % Calculate node replacement (appearance) probability.
                 logLikelihoodArr(itr) = spatialProb + angularProb;
            end
            
           % Calculate the final likelihood.
           logLikelihood = sum(logLikelihoodArr);
        end
        
%         % The function that estimates reconstruction error by comparing data features 
%         % and reconstructed features. DataType is a string that can be
%         % added to specify different cases (2D, 2.5D, 3D). 
%         function reconstructionError = CalculateReconstructionError(dataType)
%             
%         end
        
        %% Setters for the variables.
        % Sanity check function.
        function this = CheckSanity(this, arr, correctSizeArr, typeArr, nameStr)
            % Check for correct size.
            sizeArr = size(arr);
            if ~isempty(sizeArr)
                if numel(sizeArr) ~= numel(correctSizeArr)
                   error(['Error in class Estimator: ' nameStr ' is not a ' numel(size(correctSizeArr)) '-dimensional array.']);
                elseif sizeArr(2) ~= correctSizeArr(:,2)
                   error(['Error in class Estimator: The size of ' nameStr ' array is not Nx' num2str(size(correctSizeArr),2) '.']);
                end
            else
                error(['Warning in class Estimator: ' nameStr ' is empty.']);
            end
            
            % Class type check.
            if ~isa(arr, typeArr)
                error(['Error in class Estimator: ' nameStr ' is not of class ' typeArr '.']);
            end
        end
        
        % Coordinate setting function for both data and reconstructed data
        % points.
        function this = SetCoords(this, coords, typeFlag)
            % If typeFlag == 0, this function sets data points. Otherwise it
            % sets reconstructed feature points.
            if ~typeFlag
                pointString = 'dataCoords';
            else
                pointString = 'reconstructedCoords';
            end
            
            this.CheckSanity(coords, [inf, 2], 'uint16', pointString);
            this.numberOfPoints = size(coords,1);
            
            % All good, set.
           if ~typeFlag
               this.dataCoords = uint16(coords);
           else
               this.reconstructedCoords = uint16(coords);
           end
        end
        
        % Label setting function for both data and reconstructed data
        % points.
        function this = SetLabels(this, labels, typeFlag)
            % If typeFlag == 0, this function sets data labels. Otherwise it
            % sets reconstructed feature labels.
            if ~typeFlag
                labelString = 'dataLabels';
            else
                labelString = 'reconstructedLabels';
            end
            
            this.CheckSanity(labels, [inf, 1], 'uint16', labelString);
            
            % All good, set.
           if ~typeFlag
               this.dataLabels = uint16(labels);
           else
               this.reconstructedLabels = uint16(labels);
           end
        end
        
        % Depth setting function for both data and reconstructed data
        % points.
        function this = SetDepths(this, depths, typeFlag)
            % If typeFlag == 0, this function sets data depths. Otherwise it
            % sets reconstructed feature depths.
            if ~typeFlag
                depthString = 'dataDepths';
            else
                depthString = 'reconstructedDepths';
            end
            
            this.CheckSanity(depths, [inf, 1], 'uint16', depthString);
            
            % All good, set.
           if ~typeFlag
               this.dataDepths = uint16(depths);
           else
               this.reconstructedDepths = uint16(depths);
           end
        end
        
        function this = SetRadius(this, radius)
            if numel(radius) ~= 1 || ~isnumeric(radius)
                error('Error in class Estimator: radius is not a constant.');
            end
            
            this.radius = double(radius);
        end
        
        function this = SetNumberOfFeatures(this, numberOfFeatures)
            if numel(numberOfFeatures) ~= 1 || ~isnumeric(numberOfFeatures)
                error('Error in class Estimator: radius is not a constant.');
            end
            
            this.numberOfFeatures = uint8(numberOfFeatures);
        end
        
        %% Getters for the variables.
        function dataCoords = GetDataCoords(this)
           dataCoords = this.dataCoords; 
        end
        
        function dataLabels = GetDataLabels(this)
           dataLabels = this.dataLabels; 
        end
        
        function dataDepths = GetDataDepths(this)
           dataDepths = this.dataDepths;
        end
        
        function reconstructedCoords = GetReconstructedCoords(this)
           reconstructedCoords = this.reconstructedCoords; 
        end
        
        function reconstructedLabels = GetReconstructedLabels(this)
           reconstructedLabels = this.reconstructedLabels; 
        end
        
        function reconstructedDepths = GetReconstructedDepths(this)
           reconstructedDepths = this.reconstructedDepths;
        end
        
        function radius = GetRadius(this)
           radius = this.radius;
        end
        
        function radius = GetNumberOfFeatures(this)
           radius = this.numberOfFeatures;
        end
    end
    
end

