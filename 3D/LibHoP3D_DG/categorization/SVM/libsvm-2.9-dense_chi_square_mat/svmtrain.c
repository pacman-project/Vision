#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "svm.h"

#include "mex.h"
#include "svm_model_matlab.h"

#if MX_API_VER < 0x07030000
typedef int mwIndex;
#endif

#define CMD_LEN 2048
#define Malloc(type,n) (type *)malloc((n)*sizeof(type))

void print_null(const char *s) {}

void exit_with_help()
{
#ifdef _ALL_FLOAT
        mexPrintf(
        "Usage: model = svmtrain_chi2_float(training_label_vector, training_instance_matrix, 'libsvm_options');\n"
#else
	mexPrintf(
	"Usage: model = svmtrain_chi2(training_label_vector, training_instance_matrix, 'libsvm_options');\n"
#endif
	"libsvm_options:\n"
	"-s svm_type : set type of SVM (default 0)\n"
	"	0 -- C-SVC\n"
	"	1 -- nu-SVC\n"
	"	2 -- one-class SVM\n"
	"	3 -- epsilon-SVR\n"
	"	4 -- nu-SVR\n"
	"-t kernel_type : set type of kernel function (default 2)\n"
	"	0 -- linear: u'*v\n"
	"	1 -- polynomial: (gamma*u'*v + coef0)^degree\n"
	"	2 -- radial basis function: exp(-gamma*|u-v|^2)\n"
	"	3 -- sigmoid: tanh(gamma*u'*v + coef0)\n"
	"	4 -- precomputed kernel (kernel values in training_instance_matrix)\n"
	"	5 -- exponential chi-square kernel: exp(- gamma * chi^2(x,y))\n"
	"	6 -- multiple exponential chi-square kernels. sum(weight_i exp( - gamma_i * chi^2(x_i,y_i))) \n"
	"See the documentation for -u for instructions on how to specify the parameters weight_i, gamma_i, etc.\n"
	"-d degree : set degree in kernel function (default 3)\n"
	"-g gamma : set gamma in kernel function (default 1/k)\n"
	"-r coef0 : set coef0 in kernel function (default 0)\n"
	"-c cost : set the parameter C of C-SVC, epsilon-SVR, and nu-SVR (default 1)\n"
	"-n nu : set the parameter nu of nu-SVC, one-class SVM, and nu-SVR (default 0.5)\n"
	"-p epsilon : set the epsilon in loss function of epsilon-SVR (default 0.1)\n"
	"-m cachesize : set cache memory size in MB (default 100)\n"
	"-e epsilon : set tolerance of termination criterion (default 0.001)\n"
	"-h shrinking : whether to use the shrinking heuristics, 0 or 1 (default 1)\n"
	"-b probability_estimates : whether to train a SVC or SVR model for probability estimates, 0 or 1 (default 0)\n"
	"-wi weight : set the parameter C of class i to weight*C, for C-SVC (default 1)\n"
	"-v n : n-fold cross validation mode\n"
	"-q : quiet mode (no outputs)\n"
	"-u : must be put at the end of all the parameters. Only used for -t 6, multiple exponential chi-square kernels\n"
	"     Use -u length -w weight -g gamma for each kernel.\n"
        "     Example: -u 850 -w 1.2 -g 1.5 420 -w 1.3 -g 1.1 specifies 2 kernels, \n"
	"              the first one is on the first 850 dimensions of the data,\n"
        "              with a weight 1.2 and gamma 1.5, the second one is on the next 420 dimensions,\n"
	"              with a weight 1.3 and gamma 1.1\n"
	);
}

/* svm arguments*/
struct svm_parameter param;		/* set by parse_command_line */
struct svm_problem prob;		/* set by read_problem */
struct svm_model *model;
struct svm_node *x_space;
int cross_validation;
int nr_fold;

void (*svm_default_print_string) (const char *) = NULL;

SVM_REAL do_cross_validation()
{
	int i;
	int total_correct = 0;
	SVM_REAL total_error = 0;
	SVM_REAL sumv = 0, sumy = 0, sumvv = 0, sumyy = 0, sumvy = 0;
	SVM_REAL *target = Malloc(SVM_REAL,prob.l);
	SVM_REAL retval = 0.0;

	svm_cross_validation(&prob,&param,nr_fold,target);
	if(param.svm_type == EPSILON_SVR ||
	   param.svm_type == NU_SVR)
	{
		for(i=0;i<prob.l;i++)
		{
			SVM_REAL y = prob.y[i];
			SVM_REAL v = target[i];
			total_error += (v-y)*(v-y);
			sumv += v;
			sumy += y;
			sumvv += v*v;
			sumyy += y*y;
			sumvy += v*y;
		}
		mexPrintf("Cross Validation Mean squared error = %g\n",total_error/prob.l);
		mexPrintf("Cross Validation Squared correlation coefficient = %g\n",
			((prob.l*sumvy-sumv*sumy)*(prob.l*sumvy-sumv*sumy))/
			((prob.l*sumvv-sumv*sumv)*(prob.l*sumyy-sumy*sumy))
			);
		retval = total_error/prob.l;
	}
	else
	{
		for(i=0;i<prob.l;i++)
			if(target[i] == prob.y[i])
				++total_correct;
		mexPrintf("Cross Validation Accuracy = %g%%\n",100.0*total_correct/prob.l);
		retval = 100.0*total_correct/prob.l;
	}
	free(target);
	return retval;
}

/* nrhs should be 3 */
int parse_command_line(int nrhs, const mxArray *prhs[], char *model_file_name)
{
	int i, argc = 1;
	char cmd[CMD_LEN];
	char *argv[CMD_LEN/2];
#ifdef _USE_CHI_SQUARE
	const char blank[2] = " ";
#endif

	/* default values */
	param.svm_type = C_SVC;
	param.kernel_type = RBF;
	param.degree = 3;
	param.gamma = 0;	
	param.coef0 = 0;
	param.nu = 0.5;
	param.cache_size = 100;
	param.C = 1;
	param.eps = 1e-3;
	param.p = 0.1;
	param.shrinking = 1;
	param.probability = 0;
	param.nr_weight = 0;
	param.weight_label = NULL;
	param.weight = NULL;
#ifdef _USE_CHI_SQUARE
	param.multiple_kernel_parm = NULL;
	param.feat_parm_text[0] = '\0';
#endif
	cross_validation = 0;
	
	if (svm_default_print_string == NULL)
		svm_default_print_string = svm_print_string;
	else
		svm_print_string = svm_default_print_string;


	if(nrhs <= 1)
		return 1;

	if(nrhs > 2)
	{
		
		mxGetString(prhs[2], cmd, mxGetN(prhs[2]) + 1);
		if((argv[argc] = strtok(cmd, " ")) != NULL)
			while((argv[++argc] = strtok(NULL, " ")) != NULL)
				;
	}


	for(i=1;i<argc;i++)
	{
		if(argv[i][0] != '-') break;
		++i;
		if(i>=argc && argv[i-1][1] != 'q')
			return 1;
		switch(argv[i-1][1])
		{
			case 's':
				param.svm_type = atoi(argv[i]);
				break;
			case 't':
				param.kernel_type = atoi(argv[i]);
				break;
			case 'd':
				param.degree = atoi(argv[i]);
				break;
			case 'g':
				param.gamma = atof(argv[i]);
				break;
			case 'r':
				param.coef0 = atof(argv[i]);
				break;
			case 'n':
				param.nu = atof(argv[i]);
				break;
			case 'm':
				param.cache_size = atof(argv[i]);
				break;
			case 'c':
				param.C = atof(argv[i]);
				break;
			case 'e':
				param.eps = atof(argv[i]);
				break;
			case 'p':
				param.p = atof(argv[i]);
				break;
			case 'h':
				param.shrinking = atoi(argv[i]);
				break;
			case 'b':
				param.probability = atoi(argv[i]);
				break;
			case 'q':
				svm_print_string = &print_null;
				i--;
				break;
			case 'v':
				cross_validation = 1;
				nr_fold = atoi(argv[i]);
				if(nr_fold < 2)
				{
					mexPrintf("n-fold cross validation: n must >= 2\n");
					return 1;
				}
				break;
			case 'w':
				++param.nr_weight;
				param.weight_label = (int *)realloc(param.weight_label,sizeof(int)*param.nr_weight);
				param.weight = (SVM_REAL *)realloc(param.weight,sizeof(SVM_REAL)*param.nr_weight);
				param.weight_label[param.nr_weight-1] = atoi(&argv[i-1][2]);
				param.weight[param.nr_weight-1] = atof(argv[i]);
				break;
#ifdef _USE_CHI_SQUARE
			case 'u':
// Make sure that -u is the last one. Need to strcat a lot of argvs back
				new_featparm(param);
				for(;i<argc;i++)
				{
					strcat(param.feat_parm_text, argv[i]);
					strcat(param.feat_parm_text, blank);
				}
				parse_custom_arguments(param.feat_parm_text, param.multiple_kernel_parm);
				return 0;
#endif
			default:
				mexPrintf("Unknown option -%c\n", argv[i-1][1]);
				return 1;
		}
	}
	return 0;
}


int read_problem_dense(const mxArray *label_vec, const mxArray *instance_mat)
{
	int i, j, k;
	int elements, max_index, sc, label_vector_row_num;
	SVM_REAL *samples, *labels;

	prob.x = NULL;
	prob.y = NULL;
	x_space = NULL;

#ifdef _ALL_FLOAT
	labels = (float *) mxGetData(label_vec);
	samples = (float *) mxGetData(instance_mat);
#else
	labels = mxGetPr(label_vec);
	samples = mxGetPr(instance_mat);
#endif
	sc = (int)mxGetN(instance_mat);

	elements = 0;
	prob.l = (int)mxGetM(instance_mat);
	label_vector_row_num = (int)mxGetM(label_vec);

	if(label_vector_row_num!=prob.l)
	{
		mexPrintf("Length of label vector does not match # of instances.\n");
		return -1;
	}
/* No need of elements in dense format */
#ifndef _DENSE_REP
	if(param.kernel_type == PRECOMPUTED)
		elements = prob.l * (sc + 1);
	else
	{
		for(i = 0; i < prob.l; i++)
		{
			for(k = 0; k < sc; k++)
				if(samples[k * prob.l + i] != 0)
					elements++;
			
			elements++;
		}
	}
#endif

	prob.y = Malloc(SVM_REAL,prob.l);
#ifdef _DENSE_REP
	prob.x = Malloc(struct svm_node, prob.l);
#else
	prob.x = Malloc(struct svm_node *,prob.l);
	x_space = Malloc(struct svm_node, elements);
#endif

	max_index = sc;
	j = 0;
	for(i = 0; i < prob.l; i++)
	{
		prob.y[i] = labels[i];
#ifdef _DENSE_REP
		(prob.x+i)->values = Malloc(SVM_REAL,sc);
		(prob.x+i)->dim = sc;
		for(k = 0; k < sc; k++)
			(prob.x+i)->values[k] = samples[k * prob.l + i];
#else
		prob.x[i] = &x_space[j];

		for(k = 0; k < sc; k++)
		{
			if(param.kernel_type == PRECOMPUTED || samples[k * prob.l + i] != 0)
			{
				x_space[j].index = k + 1;
				x_space[j].value = samples[k * prob.l + i];
				j++;
			}
		}
		x_space[j++].index = -1;
#endif
	}

	if(param.gamma == 0 && max_index > 0)
		param.gamma = 1.0/max_index;

	if(param.kernel_type == PRECOMPUTED)
		for(i=0;i<prob.l;i++)
		{
#ifdef _DENSE_REP
			if ((int)(prob.x+i)->values[0] < 0 || (int)(prob.x+i)->values[0] > max_index)
			{
				mexPrintf("Wrong input format: sample_serial_number out of range\n");
				return -1;
			}
#else
			if ((int)prob.x[i][0].value <= 0 || (int)prob.x[i][0].value > max_index)
			{
				mexPrintf("Wrong input format: sample_serial_number out of range\n");
				return -1;
			}
#endif
		}

	return 0;
}

int read_problem_sparse(const mxArray *label_vec, const mxArray *instance_mat)
{
#ifdef _DENSE_REP
	mexPrintf("LIBSVM_Dense currently doesn't support sparse matrix input!\n");
	return -1;
#else
	int i, j, k, low, high;
	mwIndex *ir, *jc;
	int elements, max_index, num_samples, label_vector_row_num;
	SVM_REAL *samples, *labels;
	mxArray *instance_mat_col; 

	prob.x = NULL;
	prob.y = NULL;
	x_space = NULL;

	{
		mxArray *prhs[1], *plhs[1];
		prhs[0] = mxDuplicateArray(instance_mat);
		if(mexCallMATLAB(1, plhs, 1, prhs, "transpose"))
		{
			mexPrintf("Error: cannot transpose training instance matrix\n");
			return -1;
		}
		instance_mat_col = plhs[0];
		mxDestroyArray(prhs[0]);
	}

	labels = mxGetPr(label_vec);
	samples = mxGetPr(instance_mat_col);
	ir = mxGetIr(instance_mat_col);
	jc = mxGetJc(instance_mat_col);

	num_samples = (int)mxGetNzmax(instance_mat_col);

	prob.l = (int)mxGetN(instance_mat_col);
	label_vector_row_num = (int)mxGetM(label_vec);

	if(label_vector_row_num!=prob.l)
	{
		mexPrintf("Length of label vector does not match # of instances.\n");
		return -1;
	}

	elements = num_samples + prob.l;
	max_index = (int)mxGetM(instance_mat_col);

	prob.y = Malloc(SVM_REAL,prob.l);
	prob.x = Malloc(struct svm_node *,prob.l);
	x_space = Malloc(struct svm_node, elements);

	j = 0;
	for(i=0;i<prob.l;i++)
	{
		prob.x[i] = &x_space[j];
		prob.y[i] = labels[i];
		low = (int)jc[i], high = (int)jc[i+1];
		for(k=low;k<high;k++)
		{
			x_space[j].index = (int)ir[k] + 1;
			x_space[j].value = samples[k];
			j++;
	 	}
		x_space[j++].index = -1;
	}

	if(param.gamma == 0 && max_index > 0)
		param.gamma = 1.0/max_index;

	return 0;
#endif
}

static void fake_answer(mxArray *plhs[])
{
	plhs[0] = mxCreateDoubleMatrix(0, 0, mxREAL);
}

void mexFunction( int nlhs, mxArray *plhs[],
		int nrhs, const mxArray *prhs[] )
{
	const char *error_msg;
	int i;

	srand(1);

	if(nrhs > 0 && nrhs < 4)
	{
		int err;

#ifdef _ALL_FLOAT
		if(!mxIsSingle(prhs[0]) || !mxIsSingle(prhs[1])) {
			mexPrintf("Error: label vector and instance matrix must be single\n");
			fake_answer(plhs);
			return;
		}
#else
		if(!mxIsDouble(prhs[0]) || !mxIsDouble(prhs[1])) {
			mexPrintf("Error: label vector and instance matrix must be double\n");
			fake_answer(plhs);
			return;
		}
#endif

		if(parse_command_line(nrhs, prhs, NULL))
		{
			exit_with_help();
			svm_destroy_param(&param);
			fake_answer(plhs);
			return;
		}

		if(mxIsSparse(prhs[1]))
		{
			if(param.kernel_type == PRECOMPUTED)
			{
				mxArray *rhs[1], *lhs[1];

				rhs[0] = mxDuplicateArray(prhs[1]);
				if(mexCallMATLAB(1, lhs, 1, rhs, "full"))
				{
					mexPrintf("Error: cannot generate a full training instance matrix\n");
					svm_destroy_param(&param);
					fake_answer(plhs);
					return;
				}
				err = read_problem_dense(prhs[0], lhs[0]);
				mxDestroyArray(lhs[0]);
				mxDestroyArray(rhs[0]);
			}
			else
				err = read_problem_sparse(prhs[0], prhs[1]);
		}
		else
			err = read_problem_dense(prhs[0], prhs[1]);

		error_msg = svm_check_parameter(&prob, &param);

		if(err || error_msg)
		{
			if (error_msg != NULL)
				mexPrintf("Error: %s\n", error_msg);
			svm_destroy_param(&param);
			free(prob.y);
#ifdef _DENSE_REP
			for(i=0;i<prob.l;i++)
				free((prob.x+i)->values);
#endif
			free(prob.x);
			free(x_space);
			fake_answer(plhs);
			return;
		}

		if(cross_validation)
		{
			SVM_REAL *ptr;
#ifdef _ALL_FLOAT
			plhs[0] = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
			ptr = (float *) mxGetData(plhs[0]);
#else
			plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
			ptr = mxGetPr(plhs[0]);
#endif
			ptr[0] = do_cross_validation();
		}
		else
		{
			int nr_feat = (int)mxGetN(prhs[1]);
			const char *error_msg;
			model = svm_train(&prob, &param);
			error_msg = model_to_matlab_structure(plhs, nr_feat, model);
			if(error_msg)
				mexPrintf("Error: can't convert libsvm model to matrix structure: %s\n", error_msg);
			svm_destroy_model(model);
			model = 0;
		}
		svm_destroy_param(&param);
		free(prob.y);
#ifdef _DENSE_REP
		for(i=0;i<prob.l;i++)
			free((prob.x+i)->values);
#endif
		free(prob.x);
		free(x_space);
		prob.x = 0;
		prob.y = 0;
		x_space = 0;
	}
	else
	{
		exit_with_help();
		fake_answer(plhs);
		return;
	}
}
