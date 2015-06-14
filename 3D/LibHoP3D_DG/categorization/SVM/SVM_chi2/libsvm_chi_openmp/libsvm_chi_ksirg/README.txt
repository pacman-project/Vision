
Modified LibSVM with CHI^2 Kernel and exponential CHI^2, 
modifications done by Krzysztof Sopy³a (http://wmii.uwm.edu.pl/~ksopyla/libsvm_chi2)
This library contains two definition of Chi Squared kernels, which are the moste popular in publication

1) K(x,y)= 1-2*sum( (xi-yi)^2/(xi+yi)); # normal CHI^2
2) K(x,y)= sum( (xi*yi)/(xi+yi)); # normalized CHI^2
3) K(x,y)= exp(-gamma*Chi-Squared(x,y));

These kernels are avaliable through matlab interface.


Build
======================
There are already compiled .mex files for windows x64 architecture, but if you want build for different platform just write

matlab>make


Usage
======================

Two new parameter options was added, now parameter "t" can accept values equals 5,6,7

matlab> model = svmtrain(training_label_vector, training_instance_matrix, 'libsvm_options');
options:
	-s svm_type : set type of SVM (default 0)"
		0 -- C-SVC"
		1 -- nu-SVC
		2 -- one-class SVM
		3 -- epsilon-SVR
		4 -- nu-SVR
	-t kernel_type : set type of kernel function (default 2)
		0 -- linear: u'*v
		1 -- polynomial: (gamma*u'*v + coef0)^degree
		2 -- radial basis function: exp(-gamma*|u-v|^2)
		3 -- sigmoid: tanh(gamma*u'*v + coef0)
		4 -- precomputed kernel (kernel values in training_set_file)
		5 -- chi-squaree kernel k(x,y)=1-2*sum( (xi-yi)^2/(xi+yi))
		6 -- norm chi-squaree kernel k(x,y)=sum( xi*yi/(xi+yi
		7 -- exponential chi-squaree kernel k(x,y)=exp(-gamma*sum( (xi-yi)^2/(xi+yi)))
	-d degree : set degree in kernel function (default 3)
	-g gamma : set gamma in kernel function (default 1/num_features)
	-r coef0 : set coef0 in kernel function (default 0)
	-c cost : set the parameter C of C-SVC, epsilon-SVR, and nu-SVR (default 1)
	-n nu : set the parameter nu of nu-SVC, one-class SVM, and nu-SVR (default 0.5)
	-p epsilon : set the epsilon in loss function of epsilon-SVR (default 0.1)
	-m cachesize : set cache memory size in MB (default 100)
	-e epsilon : set tolerance of termination criterion (default 0.001)
	-h shrinking : whether to use the shrinking heuristics, 0 or 1 (default 1)
	-b probability_estimates : whether to train a SVC or SVR model for probability estimates, 0 or 1 (default 0)
	-wi weight : set the parameter C of class i to weight*C, for C-SVC (default 1)
	-v n: n-fold cross validation mode
	-q : quiet mode (no outputs)

Example
======================
See test_libsvm.m

Classification with chi^2 kernel, with C=4

tr_path='a1a.train';
tst_path='a1a.train';

[trYY trXX]=libsvmread(tr_path);
[tstYY tstXX]=libsvmread(tst_path);
trXXn=trXX;
%l1 - norm
trXXn=bsxfun(@rdivide,trXXn,sum(trXXn,2)); 
tic;
model = svmtrain(trYY, trXXn,'-c 4 -t 5');
modelTime=toc;

tstXXn=tstXX;
tstXXn=bsxfun(@rdivide,tstXXn,sum(tstXXn,2)); %l1 - norm

tic
[pred, acc, dec_vals] = svmpredict(tstYY, tstXXn, model);
predTime = toc;
ss=sprintf('libsvm chi^2 acc=%0.5g modeltime=%g predtime=%g \n',acc(1),modelTime, predTime);
disp(ss);


Additional Information
======================

LibSVM matlab interface was initially written by Jun-Cheng Chen, Kuan-Jen Peng,
Chih-Yuan Yang and Chih-Huai Cheng from Department of Computer
Science, National Taiwan University. The current version was prepared
by Rong-En Fan and Ting-Fan Wu. If you find this tool useful, please
cite LIBSVM as follows

Chih-Chung Chang and Chih-Jen Lin, LIBSVM : a library for support
vector machines. ACM Transactions on Intelligent Systems and
Technology, 2:27:1--27:27, 2011. Software available at
http://www.csie.ntu.edu.tw/~cjlin/libsvm

For any question, please contact Chih-Jen Lin <cjlin@csie.ntu.edu.tw>,
or check the FAQ page:

http://www.csie.ntu.edu.tw/~cjlin/libsvm/faq.html#/Q9:_MATLAB_interface
