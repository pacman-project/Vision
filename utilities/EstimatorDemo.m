% Create mock data.
dataPoints = uint16([2, 15; ...
                      10, 25; ...
                      8, 23; ...
                      15, 18; ...
                      10, 10]);
 reconstructedPoints = uint16([5, 13; ...
                                     12, 22; ...
                                     8, 20; ...
                                     13, 19; ...
                                     0, 0]);     % This row is fake data.
 dataLabels = uint16([1;1;1;2;4]);
 reconstructedLabels = uint16([1;1;2;4;3]);
 numberOfFeatures = 8;
 
 % Create object and set values.
 estimator = Estimator;
 estimator = estimator.SetCoords(dataPoints, 0);
 tmpArr = estimator.GetDataCoords();
 if ~isequal(tmpArr, dataPoints)
      error('Coordinate assignment error.');
 end
      
      
 estimator = estimator.SetCoords(reconstructedPoints,1);
 tmpArr = estimator.GetReconstructedCoords();
 if ~isequal(tmpArr, reconstructedPoints)
      error('Coordinate assignment error.');
 end
 
 estimator = estimator.SetLabels(dataLabels,0);
 tmpArr = estimator.GetDataLabels();
 if ~isequal(tmpArr, dataLabels)
      error('Label assignment error.');
 end
 
 estimator = estimator.SetLabels(reconstructedLabels,1);
 tmpArr = estimator.GetReconstructedLabels();
 if ~isequal(tmpArr, reconstructedLabels)
      error('Label assignment error.');
 end
 
 estimator = estimator.SetNumberOfFeatures(numberOfFeatures);
 tmp = estimator.GetNumberOfFeatures();
 if ~isequal(tmp, numberOfFeatures)
      error('Number of features assignment error.');
 end
 
 % Finally, calculate the probability.
 prob = estimator.CalculateDataLikelihood();
       