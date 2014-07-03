function [svmstruct] = one_vs_all_regression_libsvm_multiple_chi2(quality, data, C, par)
% One-against-all regression training using multiple chi-square kernel.
% quality: The quality matrix (e.g. overlap) on all the n examples for all the c classes (n * c array)
% data: A 1*k cell array, each cell contains the features of one kernel, with a total of k kernels
%       In each cell, the data should be an n * d matrix (d is the dimensionality of the features of this kernel)
% C  :  the SVM regularization parameter C
% par:  parameters for customizing the SVM.
%       par.cache_size: The size you are willing to spend for SVM caching on the kernel matrix 
%                        , in megabytes, e.g. cache_size = 4000 specifies a 4GB cache
%       par.weight: weight vector, 1 weight for each kernel
%       par.gamma: gamma vector (kernel width), 1 gamma for each kernel
%       par.custom_svm_option: other svm options you would like to use

    svmstruct = cell(size(quality,2),1);

    if isa(data{1},'single')
	fhandle = @svmtrain_chi2_float;
	binary_label_sets = single(quality);
    else
        fhandle = @svmtrain_chi2;
        binary_label_sets = double(quality);
    end
    
% Use the 6-th kernel (multiple chi-square) and the specified C
    options_string = ['-t 6 -c ' sprintf('%f', C)];
% Append cache size
    if isfield(par,'cache_size') && ~isempty(par.cache_size)
        options_string = [options_string '-m ' int2str(par.cache_size)];
    end
% Append -s 3 for regression
    options_string = [options_string ' -s 3'];
% Append customized options
    if isfield(par, 'custom_svm_option')
        options_string = [options_string ' ' par.custom_svm_option];
    end
% Construct the options for multiple chi-square kernels, from data and par.weight, par.gamma
    option_string2 = construct_multiple_chi_option(data, par);
    options_string = [options_string option_string2];
% Convert data into matrix form
    if size(data,1) > size(data,2)
        data = data';
    end
    tdata = cell2mat(data);
    disp(['LIBSVM options: ' options_string]);
% Train SVM one by one
    for i=1:size(binary_label_sets, 2)
        t =tic();
        svmstruct{i} = fhandle(binary_label_sets(:,i), tdata, options_string);
        toc(t);
    end
end

% Note that the chi-square options must be at the end of the option string,
% because it has multiple options after -u
function options_string = construct_multiple_chi_option(data, par)

    options_string = ' -u';
    for i=1:length(data)
        options_string = [options_string sprintf(' %d -w %.4f -g %.4f', size(data{i},2), par.weights(i), par.gamma(i))];
    end
end
