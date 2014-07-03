function [predicted_labels, pred] = one_vs_all_predict_libsvm_multiple_chi2(test_qualities, data, model)
% One-against-all regression prediction.
% test_qualities: The quality matrix (e.g. overlap) on all the n examples for all the c classes (n * c array)
% data: A 1*k cell array, each cell contains the features of one kernel, with a total of k kernels
%       In each cell, the data should be an n * d matrix (d is the dimensionality of the features of this kernel)
% model: The model trained by one_vs_all_train_libsvm_multiple_chi2
    pred = cell(length(model),1);
    tdata = cell2mat(data);
    real_size = [size(tdata,1) length(model)];
    if isempty(test_qualities) || sum(size(test_qualities) ~= real_size)
        test_qualities = zeros(size(tdata,1), length(model));
    end
    % Usage depends on whether the model is single or double
	if isa(model{1}.sv_coef, 'single')
    	fhandle = @svmpredict_chi2_float;
        tdata = single(tdata);
        test_qualities = single(test_qualities);
    else
		fhandle = @svmpredict_chi2;
        tdata = double(tdata);
        test_qualities = double(test_qualities);
	end
    % No parfor here because parallelization is done inside svmpredict_chi2
    % already.
    for i=1:length(model)
        t = tic();
        [~,~, pred{i}] = fhandle(test_qualities(:,i), tdata, model{i});
        if(~isempty(model{i}.Label) && model{i}.Label(1) == -1)
            pred{i} = -pred{i};
        end
        toc(t);
    end
    
    pred = cell2mat(pred');
    [~, predicted_labels] = max(pred, [], 2);
end