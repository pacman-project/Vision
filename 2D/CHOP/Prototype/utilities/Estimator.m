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
        dataCoords@uint16; % Nx2 array of uint16. Each row holds coordinates (x,y) of a single data point.
        dataLabels@uint16; % Nx1 array of uint16. Point labels.
        dataDepths@uint16; % Nx1 array of uint16. Used for 2.5D. Point depths.
        reconstructedCoords@uint16; % Nx2 array of uint16. Each row holds coordinates (x,y) of a reconstructed feature.
        reconstructedLabels@uint16; % Nx1 array of uint16. Reconstructed feature labels.
        reconstructedDepths@uint16; % Nx1 array of uint16. Used for 2.5D. Reconstructed feature depths.
        radius@double;      % Constant. The radius in which the algorithm looks for correspondances.
    end
    
    methods
        %% Core methods
        % The function that estimates data likelihood by comparing data features 
        % and reconstructed features. DataType is a string that can be
        % added to specify different cases (2D, 2.5D, 3D). 
        function likelihood = CalculateDataLikelihood(dataType)
            
        end
        
        
        % The function that estimates reconstruction error by comparing data features 
        % and reconstructed features. DataType is a string that can be
        % added to specify different cases (2D, 2.5D, 3D). 
        function reconstructionError = CalculateReconstructionError(dataType)
            
        end
        
        %% Setters for the variables.
        % Sanity check function.
        function isSane = CheckSanity(arr, correctSizeArr, typeArr, nameStr)
            % Check for correct size.
            isSane = true;
            sizeArr = size(arr);
            if ~isempty(sizeArr)
                if numel(sizeArr) ~= numel(correctSizeArr)
                   isSane = false;
                   error(['Error in class Estimator: ' nameStr ' is not a ' numel(correctSizeArr) '-dimensional array.']);
                elseif sizeArr(2) ~= 2
                   isSane = false;
                   error(['Error in class Estimator: The size of ' nameStr ' array is not Nx' num2str(correctSizeArr,2) '.']);
                end
            else
                warning(['Warning in class Estimator: ' nameStr ' is empty.']);
            end
            
            % Class type check.
            if ~isa(coords, typeArr)
                isSane = false;
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
            
            isSane = CheckSanity(coords, [inf, 2], 'uint16', pointString);
            
            % All good, set.
            if isSane
                if ~typeFlag
                    this.dataCoords = coords;
                else
                    this.reconstructedCoords = coords;
                end
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
            
            isSane = CheckSanity(labels, [inf, 1], 'uint16', labelString);
            
            % All good, set.
            if isSane
                if ~typeFlag
                    this.dataLabels = labels;
                else
                    this.reconstructedLabels = labels;
                end
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
            
            isSane = CheckSanity(depths, [inf, 1], 'uint16', depthString);
            
            % All good, set.
            if isSane
                if ~typeFlag
                    this.dataDepths = depths;
                else
                    this.reconstructedDepths = depths;
                end
            end
        end
        
        function this = SetRadius(this, radius)
            if numel(radius) ~= 1 || ~isnumeric(radius)
                error('Error in class Estimator: radius is not a constant.');
            end
            
            this.radius = double(radius);
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
    end
    
end

