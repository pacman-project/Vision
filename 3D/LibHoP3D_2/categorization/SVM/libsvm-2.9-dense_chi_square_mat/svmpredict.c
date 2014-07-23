#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "svm.h"

#include "mex.h"
#include "svm_model_matlab.h"

#if MX_API_VER < 0x07030000
typedef int mwIndex;
#endif

#define CMD_LEN 2048

void read_sparse_instance(const mxArray *prhs, int index, struct svm_node *x)
{
#ifdef _DENSE_REP
	mexPrintf("Sparse instance for LIBSVM_Dense is not supported currently!");
#else
	int i, j, low, high;
	mwIndex *ir, *jc;
	SVM_REAL *samples;

	ir = mxGetIr(prhs);
	jc = mxGetJc(prhs);
	samples = mxGetPr(prhs);

	j = 0;
	low = (int)jc[index], high = (int)jc[index+1];
	for(i=low;i<high;i++)
	{
		x[j].index = (int)ir[i] + 1;
		x[j].value = samples[i];
		j++;
	}
	x[j].index = -1;
#endif
}

static void fake_answer(mxArray *plhs[])
{
	plhs[0] = mxCreateDoubleMatrix(0, 0, mxREAL);
	plhs[1] = mxCreateDoubleMatrix(0, 0, mxREAL);
	plhs[2] = mxCreateDoubleMatrix(0, 0, mxREAL);
}

void predict(mxArray *plhs[], const mxArray *prhs[], struct svm_model *model, const int predict_probability)
{
	int label_vector_row_num, label_vector_col_num;
	int feature_number, testing_instance_number;
	int instance_index;
	SVM_REAL *ptr_instance, *ptr_label, *ptr_predict_label; 
	SVM_REAL *ptr_prob_estimates, *ptr_dec_values, *ptr;
	struct svm_node *x;
	mxArray *pplhs[1];

	int correct = 0;
	int total = 0;
	SVM_REAL error = 0;
	SVM_REAL sump = 0, sumt = 0, sumpp = 0, sumtt = 0, sumpt = 0;

	int svm_type=svm_get_svm_type(model);
	int nr_class=svm_get_nr_class(model);
	SVM_REAL *prob_estimates=NULL;

	feature_number = (int)mxGetN(prhs[1]);
	testing_instance_number = (int)mxGetM(prhs[1]);
	label_vector_row_num = (int)mxGetM(prhs[0]);
	label_vector_col_num = (int)mxGetN(prhs[0]);

	if(label_vector_row_num!=testing_instance_number)
	{
		mexPrintf("Length of label vector does not match # of instances.\n");
		fake_answer(plhs);
		return;
	}
	if(label_vector_col_num!=1)
	{
		mexPrintf("label (1st argument) should be a vector (# of column is 1).\n");
		fake_answer(plhs);
		return;
	}

#ifdef _ALL_FLOAT
	ptr_instance = (float *)mxGetData(prhs[1]);
	ptr_label    = (float *)mxGetData(prhs[0]);
#else
	ptr_instance = mxGetPr(prhs[1]);
	ptr_label    = mxGetPr(prhs[0]);
#endif

	if(mxIsSparse(prhs[1]))
	{
		if(model->param.kernel_type == PRECOMPUTED)
		{
			
			mxArray *rhs[1], *lhs[1];
			rhs[0] = mxDuplicateArray(prhs[1]);
			if(mexCallMATLAB(1, lhs, 1, rhs, "full"))
			{
				mexPrintf("Error: cannot full testing instance matrix\n");
				fake_answer(plhs);
				return;
			}
#ifdef _ALL_FLOAT
			ptr_instance = (float *)mxGetData(lhs[0]);
#else
			ptr_instance = mxGetPr(lhs[0]);
#endif
			mxDestroyArray(rhs[0]);
		}
		else
		{
			mxArray *pprhs[1];
			pprhs[0] = mxDuplicateArray(prhs[1]);
			if(mexCallMATLAB(1, pplhs, 1, pprhs, "transpose"))
			{
				mexPrintf("Error: cannot transpose testing instance matrix\n");
				fake_answer(plhs);
				return;
			}
		}
	}

	if(predict_probability)
	{
		if(svm_type==NU_SVR || svm_type==EPSILON_SVR)
			mexPrintf("Prob. model for test data: target value = predicted value + z,\nz: Laplace distribution e^(-|z|/sigma)/(2sigma),sigma=%g\n",svm_get_svr_probability(model));
		else
			prob_estimates = (SVM_REAL *) malloc(nr_class*sizeof(SVM_REAL));
	}

#ifdef _ALL_FLOAT
	plhs[0] = mxCreateNumericMatrix(testing_instance_number, 1, mxSINGLE_CLASS, mxREAL);
	if(predict_probability)
	{
		
		if(svm_type==C_SVC || svm_type==NU_SVC)
			plhs[2] = mxCreateNumericMatrix(testing_instance_number, nr_class, mxSINGLE_CLASS,mxREAL);
		else
			plhs[2] = mxCreateNumericMatrix(0, 0, mxSINGLE_CLASS,mxREAL);
	}
	else
	{
		
		if(svm_type == ONE_CLASS ||
		   svm_type == EPSILON_SVR ||
		   svm_type == NU_SVR)
			plhs[2] = mxCreateNumericMatrix(testing_instance_number, 1, mxSINGLE_CLASS,mxREAL);
		else
			plhs[2] = mxCreateNumericMatrix(testing_instance_number, nr_class*(nr_class-1)/2, mxSINGLE_CLASS,mxREAL);
	}

	ptr_predict_label = (float *)mxGetData(plhs[0]);
	ptr_prob_estimates = (float *)mxGetData(plhs[2]);
	ptr_dec_values = (float *)mxGetData(plhs[2]);
#else
	plhs[0] = mxCreateDoubleMatrix(testing_instance_number, 1, mxREAL);
	if(predict_probability)
	{
		
		if(svm_type==C_SVC || svm_type==NU_SVC)
			plhs[2] = mxCreateDoubleMatrix(testing_instance_number, nr_class, mxREAL);
		else
			plhs[2] = mxCreateDoubleMatrix(0, 0, mxREAL);
	}
	else
	{
		
		if(svm_type == ONE_CLASS ||
		   svm_type == EPSILON_SVR ||
		   svm_type == NU_SVR)
			plhs[2] = mxCreateDoubleMatrix(testing_instance_number, 1, mxREAL);
		else
			plhs[2] = mxCreateDoubleMatrix(testing_instance_number, nr_class*(nr_class-1)/2, mxREAL);
	}

	ptr_predict_label = mxGetPr(plhs[0]);
	ptr_prob_estimates = mxGetPr(plhs[2]);
	ptr_dec_values = mxGetPr(plhs[2]);
#endif
#ifdef _DENSE_REP
	x = (struct svm_node *)malloc(sizeof(struct svm_node));
	x->values = (SVM_REAL *) malloc(sizeof(SVM_REAL)*feature_number);
	x->dim = feature_number;
#else
	x = (struct svm_node*)malloc((feature_number+1)*sizeof(struct svm_node) );
#endif
	for(instance_index=0;instance_index<testing_instance_number;instance_index++)
	{
		int i;
		SVM_REAL target_label, predict_label;

		target_label = ptr_label[instance_index];

		if(mxIsSparse(prhs[1]) && model->param.kernel_type != PRECOMPUTED)
			read_sparse_instance(pplhs[0], instance_index, x);
		else
		{
#ifdef _DENSE_REP
			for(i=0;i<feature_number;i++)
				x->values[i] = ptr_instance[testing_instance_number*i+instance_index];
#else
			for(i=0;i<feature_number;i++)
			{
				x[i].index = i+1;
				x[i].value = ptr_instance[testing_instance_number*i+instance_index];
			}
			x[feature_number].index = -1;
#endif
		}
		if(predict_probability)
		{
			if(svm_type==C_SVC || svm_type==NU_SVC)
			{
				predict_label = svm_predict_probability(model, x, prob_estimates);
				ptr_predict_label[instance_index] = predict_label;
				for(i=0;i<nr_class;i++)
					ptr_prob_estimates[instance_index + i * testing_instance_number] = prob_estimates[i];
			} else {
				predict_label = svm_predict(model,x);
				ptr_predict_label[instance_index] = predict_label;
			}
		}
		else
		{
			predict_label = svm_predict(model,x);
			ptr_predict_label[instance_index] = predict_label;

			if(svm_type == ONE_CLASS ||
			   svm_type == EPSILON_SVR ||
			   svm_type == NU_SVR)
			{
				SVM_REAL res;
				svm_predict_values(model, x, &res);
				ptr_dec_values[instance_index] = res;
			}
			else
			{
				SVM_REAL *dec_values = (SVM_REAL *) malloc(sizeof(SVM_REAL) * nr_class*(nr_class-1)/2);
				svm_predict_values(model, x, dec_values);
				for(i=0;i<(nr_class*(nr_class-1))/2;i++)
					ptr_dec_values[instance_index + i * testing_instance_number] = dec_values[i];
				free(dec_values);
			}
		}

		if(predict_label == target_label)
			++correct;
		error += (predict_label-target_label)*(predict_label-target_label);
		sump += predict_label;
		sumt += target_label;
		sumpp += predict_label*predict_label;
		sumtt += target_label*target_label;
		sumpt += predict_label*target_label;
		++total;
	}
	if(svm_type==NU_SVR || svm_type==EPSILON_SVR)
	{
		mexPrintf("Mean squared error = %g (regression)\n",(SVM_REAL)error/(SVM_REAL)total);
		mexPrintf("Squared correlation coefficient = %g (regression)\n",
			(((SVM_REAL)total*sumpt-sump*sumt)*((SVM_REAL)total*sumpt-sump*sumt))/
			(((SVM_REAL)total*sumpp-sump*sump)*((SVM_REAL)total*sumtt-sumt*sumt))
			);
	}
	else
		mexPrintf("Accuracy = %g%% (%d/%d) (classification)\n",
			(SVM_REAL)correct/(SVM_REAL)total*100,correct,total);

#ifdef _ALL_FLOAT
	plhs[1] = mxCreateNumericMatrix(3, 1, mxSINGLE_CLASS,mxREAL);
	ptr = (float *) mxGetData(plhs[1]);
#else
	plhs[1] = mxCreateDoubleMatrix(3, 1, mxREAL);
	ptr = mxGetPr(plhs[1]);
#endif
	ptr[0] = (SVM_REAL)(correct/total*100);
	ptr[1] = (SVM_REAL)error / (SVM_REAL)total;
	ptr[2] = (SVM_REAL)((((SVM_REAL)total*sumpt-sump*sumt)*((SVM_REAL)total*sumpt-sump*sumt))/
				(((SVM_REAL)total*sumpp-sump*sump)*((SVM_REAL)total*sumtt-sumt*sumt)));
#ifdef _DENSE_REP
	free(x->values);
#endif
	free(x);
	if(prob_estimates != NULL)
		free(prob_estimates);
}

void exit_with_help()
{
#ifdef _ALL_FLOAT
        mexPrintf(
                "Usage: [predicted_label, accuracy, decision_values/prob_estimates] = svmpredict_chi2_float(testing_label_vector, testing_instance_matrix, model, 'libsvm_options')\n"
#else
	mexPrintf(
		"Usage: [predicted_label, accuracy, decision_values/prob_estimates] = svmpredict_chi2(testing_label_vector, testing_instance_matrix, model, 'libsvm_options')\n"
#endif
		"Parameters:\n"
		"  model: SVM model structure from svmtrain.\n"
		"  libsvm_options:\n"
		"    -b probability_estimates: whether to predict probability estimates, 0 or 1 (default 0); one-class SVM not supported yet\n"
		"Returns:\n"
		"  predicted_label: SVM prediction output vector.\n"
		"  accuracy: a vector with accuracy, mean squared error, squared correlation coefficient.\n"
		"  prob_estimates: If selected, probability estimate vector.\n"
	);
}

void mexFunction( int nlhs, mxArray *plhs[],
		 int nrhs, const mxArray *prhs[] )
{
	int prob_estimate_flag = 0;
	struct svm_model *model;

	if(nrhs > 4 || nrhs < 3)
	{
		exit_with_help();
		fake_answer(plhs);
		return;
	}

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

	if(mxIsStruct(prhs[2]))
	{
		const char *error_msg;

		
		if(nrhs==4)
		{
			int i, argc = 1;
			char cmd[CMD_LEN], *argv[CMD_LEN/2];

			
			mxGetString(prhs[3], cmd,  mxGetN(prhs[3]) + 1);
			if((argv[argc] = strtok(cmd, " ")) != NULL)
				while((argv[++argc] = strtok(NULL, " ")) != NULL)
					;

			for(i=1;i<argc;i++)
			{
				if(argv[i][0] != '-') break;
				if(++i>=argc)
				{
					exit_with_help();
					fake_answer(plhs);
					return;
				}
				switch(argv[i-1][1])
				{
					case 'b':
						prob_estimate_flag = atoi(argv[i]);
						break;
					default:
						mexPrintf("Unknown option: -%c\n", argv[i-1][1]);
						exit_with_help();
						fake_answer(plhs);
						return;
				}
			}
		}

		model = matlab_matrix_to_model(prhs[2], &error_msg);
		if (model == NULL)
		{
			mexPrintf("Error: can't read model: %s\n", error_msg);
			fake_answer(plhs);
			return;
		}

		if(prob_estimate_flag)
		{
			if(svm_check_probability_model(model)==0)
			{
				mexPrintf("Model does not support probabiliy estimates\n");
				fake_answer(plhs);
				svm_destroy_model(model);
				return;
			}
		}
		else
		{
			if(svm_check_probability_model(model)!=0)
				printf("Model supports probability estimates, but disabled in predicton.\n");
		}

		predict(plhs, prhs, model, prob_estimate_flag);
		
		svm_destroy_model(model);
	}
	else
	{
		mexPrintf("model file should be a struct array\n");
		fake_answer(plhs);
	}

	return;
}
