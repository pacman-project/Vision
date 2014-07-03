setenv('OMP_NUM_THREADS','2');
var = load('bow_features_train.mat');
data{1} = var.Feats;
quality = var.quality;
var = load('hog_features_train.mat');
data{2} = var.Feats;
par.gamma = [1.5 1.4758];
par.weights = [2 1];
model = one_vs_all_regression_libsvm_multiple_chi2(quality, data, 25, par);
var = load('bow_features_test.mat');
test_qualities = var.quality;
[~,test_labels] = max(test_qualities,[],2);
test_data{1} = var.Feats;
var = load('hog_features_test.mat');
test_data{2} = var.Feats;
[predicted_labels, pred] = one_vs_all_predict_libsvm_multiple_chi2(test_qualities, test_data, model);
disp('Prediction Accuracy: ');
sum(predicted_labels == test_labels) / length(test_labels)
