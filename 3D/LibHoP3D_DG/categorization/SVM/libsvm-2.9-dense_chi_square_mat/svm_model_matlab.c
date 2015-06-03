#include <stdlib.h>
#include <string.h>
#include "svm.h"

#include "mex.h"

#if MX_API_VER < 0x07030000
typedef int mwIndex;
#endif

#ifdef _USE_CHI_SQUARE
#define NUM_OF_RETURN_FIELD 11
#else
#define NUM_OF_RETURN_FIELD 10
#endif

#define Malloc(type,n) (type *)malloc((n)*sizeof(type))

static const char *field_names[] = {
	"Parameters",
	"nr_class",
	"totalSV",
	"rho",
	"Label",
	"ProbA",
	"ProbB",
	"nSV",
	"sv_coef",
	"SVs"
#ifdef _USE_CHI_SQUARE
	,"multiple_gamma"
#endif
};

const char *model_to_matlab_structure(mxArray *plhs[], int num_of_feature, struct svm_model *model)
{
	int i, j, n;
	SVM_REAL *ptr;
	mxArray *return_model, **rhs;
	int out_id = 0;

	rhs = (mxArray **)mxMalloc(sizeof(mxArray *)*NUM_OF_RETURN_FIELD);

	/* Parameters */
#ifdef _ALL_FLOAT
	rhs[out_id] = mxCreateNumericMatrix(5, 1, mxSINGLE_CLASS, mxREAL);
	ptr = (float *) mxGetData(rhs[out_id]);
#else
	rhs[out_id] = mxCreateDoubleMatrix(5, 1, mxREAL);
	ptr = mxGetPr(rhs[out_id]);
#endif
	ptr[0] = (SVM_REAL)model->param.svm_type;
	ptr[1] = (SVM_REAL)model->param.kernel_type;
	ptr[2] = (SVM_REAL)model->param.degree;
	ptr[3] = model->param.gamma;
	ptr[4] = model->param.coef0;
	out_id++;

	/* nr_class */
#ifdef _ALL_FLOAT
	rhs[out_id] = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS,mxREAL);
	ptr = (float *) mxGetData(rhs[out_id]);
#else
	rhs[out_id] = mxCreateDoubleMatrix(1, 1, mxREAL);
	ptr = mxGetPr(rhs[out_id]);
#endif
	ptr[0] = (SVM_REAL)model->nr_class;
	out_id++;

	/* total SV */
#ifdef _ALL_FLOAT
	rhs[out_id] = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS,mxREAL);
	ptr = (float *) mxGetData(rhs[out_id]);
#else
	rhs[out_id] = mxCreateDoubleMatrix(1, 1, mxREAL);
	ptr = mxGetPr(rhs[out_id]);
#endif
	ptr[0] = (SVM_REAL)model->l;
	out_id++;

	/* rho */
	n = model->nr_class*(model->nr_class-1)/2;
#ifdef _ALL_FLOAT
	rhs[out_id] = mxCreateNumericMatrix(n, 1, mxSINGLE_CLASS,mxREAL);
	ptr = (float *) mxGetData(rhs[out_id]);
#else
	rhs[out_id] = mxCreateDoubleMatrix(n, 1, mxREAL);
	ptr = mxGetPr(rhs[out_id]);
#endif
	for(i = 0; i < n; i++)
		ptr[i] = model->rho[i];
	out_id++;

	/* Label */
	if(model->label)
	{
#ifdef _ALL_FLOAT
		rhs[out_id] = mxCreateNumericMatrix(model->nr_class, 1, mxSINGLE_CLASS,mxREAL);
		ptr = (float *) mxGetData(rhs[out_id]);
#else
		rhs[out_id] = mxCreateDoubleMatrix(model->nr_class, 1, mxREAL);
		ptr = mxGetPr(rhs[out_id]);
#endif
		for(i = 0; i < model->nr_class; i++)
			ptr[i] = (SVM_REAL)model->label[i];
	}
	else
#ifdef _ALL_FLOAT
		rhs[out_id] = mxCreateNumericMatrix(0, 0, mxSINGLE_CLASS,mxREAL);
#else
		rhs[out_id] = mxCreateDoubleMatrix(0, 0, mxREAL);
#endif
	out_id++;

	/* probA */
	if(model->probA != NULL)
	{
#ifdef _ALL_FLOAT
		rhs[out_id] = mxCreateNumericMatrix(n, 1, mxSINGLE_CLASS,mxREAL);
		ptr = (float *) mxGetData(rhs[out_id]);
#else
		rhs[out_id] = mxCreateDoubleMatrix(n, 1, mxREAL);
		ptr = mxGetPr(rhs[out_id]);
#endif
		for(i = 0; i < n; i++)
			ptr[i] = model->probA[i];
	}
	else
#ifdef _ALL_FLOAT
		rhs[out_id] = mxCreateNumericMatrix(0, 0, mxSINGLE_CLASS,mxREAL);
#else
		rhs[out_id] = mxCreateDoubleMatrix(0, 0, mxREAL);
#endif
	out_id ++;

	/* probB */
	if(model->probB != NULL)
	{
#ifdef _ALL_FLOAT
		rhs[out_id] = mxCreateNumericMatrix(n, 1, mxSINGLE_CLASS,mxREAL);
		ptr = (float *) mxGetData(rhs[out_id]);
#else
		rhs[out_id] = mxCreateDoubleMatrix(n, 1, mxREAL);
		ptr = mxGetPr(rhs[out_id]);
#endif
		for(i = 0; i < n; i++)
			ptr[i] = model->probB[i];
	}
	else
#ifdef _ALL_FLOAT
		rhs[out_id] = mxCreateNumericMatrix(0, 0, mxSINGLE_CLASS,mxREAL);
#else
		rhs[out_id] = mxCreateDoubleMatrix(0, 0, mxREAL);
#endif
	out_id++;

	/* nSV */
	if(model->nSV)
	{
#ifdef _ALL_FLOAT
		rhs[out_id] = mxCreateNumericMatrix(model->nr_class, 1, mxSINGLE_CLASS,mxREAL);
		ptr = (float *) mxGetData(rhs[out_id]);
#else
		rhs[out_id] = mxCreateDoubleMatrix(model->nr_class, 1, mxREAL);
		ptr = mxGetPr(rhs[out_id]);
#endif
		for(i = 0; i < model->nr_class; i++)
			ptr[i] = (SVM_REAL)model->nSV[i];
	}
	else
#ifdef _ALL_FLOAT
		rhs[out_id] = mxCreateNumericMatrix(0, 0, mxSINGLE_CLASS,mxREAL);
#else
		rhs[out_id] = mxCreateDoubleMatrix(0, 0, mxREAL);
#endif
	out_id++;

	/* sv_coef */
#ifdef _ALL_FLOAT
	rhs[out_id] = mxCreateNumericMatrix(model->l, model->nr_class-1, mxSINGLE_CLASS,mxREAL);
	ptr = (float *) mxGetData(rhs[out_id]);
#else
	rhs[out_id] = mxCreateDoubleMatrix(model->l, model->nr_class-1, mxREAL);
	ptr = mxGetPr(rhs[out_id]);
#endif
	for(i = 0; i < model->nr_class-1; i++)
		for(j = 0; j < model->l; j++)
			ptr[(i*(model->l))+j] = model->sv_coef[i][j];
	out_id++;

	/* SVs */
	{
#ifdef _DENSE_REP
#ifdef _ALL_FLOAT
		rhs[out_id] = mxCreateNumericMatrix(model->l, num_of_feature, mxSINGLE_CLASS,mxREAL);
		ptr = (float *) mxGetData(rhs[out_id]);
#else
		rhs[out_id] = mxCreateDoubleMatrix(model->l, num_of_feature,  mxREAL);
		ptr = mxGetPr(rhs[out_id]);
#endif
		for(i=0;i < num_of_feature;i++)
			for(j=0;j < model->l;j++)
				ptr[(i*(model->l))+j] = model->SV[j].values[i];
#else
		int ir_index, nonzero_element;
		mwIndex *ir, *jc;
		mxArray *pprhs[1], *pplhs[1];	
		if(model->param.kernel_type == PRECOMPUTED)
		{
			nonzero_element = model->l;
			num_of_feature = 1;
		}
		else
		{
			nonzero_element = 0;
			for(i = 0; i < model->l; i++) {
				j = 0;
				while(model->SV[i][j].index != -1) 
				{
					nonzero_element++;
					j++;
				}
			}
		}

		/* SV in column, easier accessing */
		double *ptr2;
		rhs[out_id] = mxCreateSparse(num_of_feature, model->l, nonzero_element, mxREAL);
		ptr2 = mxGetPr(rhs[out_id]);
		ir = mxGetIr(rhs[out_id]);
		jc = mxGetJc(rhs[out_id]);
		jc[0] = ir_index = 0;		
		for(i = 0;i < model->l; i++)
		{
			if(model->param.kernel_type == PRECOMPUTED)
			{
				/* make a (1 x model->l) matrix */
				ir[ir_index] = 0; 
				ptr2[ir_index] = model->SV[i][0].value;
				ir_index++;
				jc[i+1] = jc[i] + 1;
			}
			else
			{
				int x_index = 0;
				while (model->SV[i][x_index].index != -1)
				{
					ir[ir_index] = model->SV[i][x_index].index - 1; 
					ptr2[ir_index] = model->SV[i][x_index].value;
					ir_index++, x_index++;
				}
				jc[i+1] = jc[i] + x_index;
			}
		}
		/* transpose back to SV in row */
		pprhs[0] = rhs[out_id];
		if(mexCallMATLAB(1, pplhs, 1, pprhs, "transpose"))
			return "cannot transpose SV matrix";
		rhs[out_id] = pplhs[0];
#endif
		out_id++;
	}
#ifdef _USE_CHI_SQUARE
	rhs[out_id] = mxCreateString(model->param.feat_parm_text);
#endif

	/* Create a struct matrix contains NUM_OF_RETURN_FIELD fields */
	return_model = mxCreateStructMatrix(1, 1, NUM_OF_RETURN_FIELD, field_names);

	/* Fill struct matrix with input arguments */
	for(i = 0; i < NUM_OF_RETURN_FIELD; i++)
		mxSetField(return_model,0,field_names[i],mxDuplicateArray(rhs[i]));
	/* return */
	plhs[0] = return_model;
	mxFree(rhs);

	return NULL;
}

struct svm_model *matlab_matrix_to_model(const mxArray *matlab_struct, const char **msg)
{
	int i, j, n, num_of_fields;
	SVM_REAL *ptr;
	int id = 0;
	struct svm_node *x_space;
	struct svm_model *model;
	mxArray **rhs;

	num_of_fields = mxGetNumberOfFields(matlab_struct);
	if(num_of_fields != NUM_OF_RETURN_FIELD) 
	{
		*msg = "number of return field is not correct";
		return NULL;
	}
	rhs = (mxArray **) mxMalloc(sizeof(mxArray *)*num_of_fields);

	for(i=0;i<num_of_fields;i++)
		rhs[i] = mxGetFieldByNumber(matlab_struct, 0, i);

	model = Malloc(struct svm_model, 1);
	model->rho = NULL;
	model->probA = NULL;
	model->probB = NULL;
	model->label = NULL;
	model->nSV = NULL;
	model->free_sv = 1; 

#ifdef _USE_CHI_SQUARE
	model->param.multiple_kernel_parm = NULL;
	model->param.feat_parm_text[0] = NULL;
#endif

#ifdef _ALL_FLOAT
	ptr = (float *)mxGetData(rhs[id]);
#else
	ptr = mxGetPr(rhs[id]);
#endif
	model->param.svm_type = (int)ptr[0];
	model->param.kernel_type  = (int)ptr[1];
	model->param.degree	  = (int)ptr[2];
	model->param.gamma	  = ptr[3];
	model->param.coef0	  = ptr[4];
	id++;

#ifdef _ALL_FLOAT
	ptr = (float *)mxGetData(rhs[id]);
#else
	ptr = mxGetPr(rhs[id]);
#endif
	model->nr_class = (int)ptr[0];
	id++;

#ifdef _ALL_FLOAT
	ptr = (float *)mxGetData(rhs[id]);
#else
	ptr = mxGetPr(rhs[id]);
#endif
	model->l = (int)ptr[0];
	id++;

	/* rho */
	n = model->nr_class * (model->nr_class-1)/2;
	model->rho = (SVM_REAL*) malloc(n*sizeof(SVM_REAL));
#ifdef _ALL_FLOAT
	ptr = (float *)mxGetData(rhs[id]);
#else
	ptr = mxGetPr(rhs[id]);
#endif
	for(i=0;i<n;i++)
		model->rho[i] = ptr[i];
	id++;

	/* label */
	if(mxIsEmpty(rhs[id]) == 0)
	{
		model->label = (int*) malloc(model->nr_class*sizeof(int));
#ifdef _ALL_FLOAT
	ptr = (float *)mxGetData(rhs[id]);
#else
	ptr = mxGetPr(rhs[id]);
#endif
		for(i=0;i<model->nr_class;i++)
			model->label[i] = (int)ptr[i];
	}
	id++;

	/* probA */
	if(mxIsEmpty(rhs[id]) == 0)
	{
		model->probA = (SVM_REAL*) malloc(n*sizeof(SVM_REAL));
#ifdef _ALL_FLOAT
	ptr = (float *)mxGetData(rhs[id]);
#else
	ptr = mxGetPr(rhs[id]);
#endif
		for(i=0;i<n;i++)
			model->probA[i] = ptr[i];
	}
	id++;

	/* probB */
	if(mxIsEmpty(rhs[id]) == 0)
	{
		model->probB = (SVM_REAL*) malloc(n*sizeof(SVM_REAL));
#ifdef _ALL_FLOAT
	ptr = (float *)mxGetData(rhs[id]);
#else
	ptr = mxGetPr(rhs[id]);
#endif
		for(i=0;i<n;i++)
			model->probB[i] = ptr[i];
	}
	id++;

	/* nSV */
	if(mxIsEmpty(rhs[id]) == 0)
	{
		model->nSV = (int*) malloc(model->nr_class*sizeof(int));
#ifdef _ALL_FLOAT
	ptr = (float *)mxGetData(rhs[id]);
#else
	ptr = mxGetPr(rhs[id]);
#endif
		for(i=0;i<model->nr_class;i++)
			model->nSV[i] = (int)ptr[i];
	}
	id++;

	/* sv_coef */
#ifdef _ALL_FLOAT
	ptr = (float *)mxGetData(rhs[id]);
#else
	ptr = mxGetPr(rhs[id]);
#endif
	model->sv_coef = (SVM_REAL**) malloc((model->nr_class-1)*sizeof(SVM_REAL *));
	for( i=0 ; i< model->nr_class -1 ; i++ )
		model->sv_coef[i] = (SVM_REAL*) malloc((model->l)*sizeof(SVM_REAL));
	for(i = 0; i < model->nr_class - 1; i++)
		for(j = 0; j < model->l; j++)
			model->sv_coef[i][j] = ptr[i*(model->l)+j];
	id++;

	/* SV */
	{
#ifdef _DENSE_REP
		int sr, sc;
#ifdef _ALL_FLOAT
	ptr = (float *)mxGetData(rhs[id]);
#else
	ptr = mxGetPr(rhs[id]);
#endif
		sr = (int)mxGetM(rhs[id]);
		sc = (int)mxGetN(rhs[id]);
		model->SV = (struct svm_node *) malloc(sr * sizeof(struct svm_node));
		x_space = NULL;
		for(j=0;j < model->l;j++)
		{
			model->SV[j].values = (SVM_REAL *)malloc(sc * sizeof(SVM_REAL));
			model->SV[j].dim = sc;
			for(i=0;i < sc;i++)
				model->SV[j].values[i] = ptr[(i*(model->l))+j];
		}
#else
		int sr, sc, elements;
		int num_samples;
		mwIndex *ir, *jc;
		mxArray *pprhs[1], *pplhs[1];

		/* transpose SV */
		pprhs[0] = rhs[id];
		if(mexCallMATLAB(1, pplhs, 1, pprhs, "transpose")) 
		{
			svm_destroy_model(model);
			*msg = "cannot transpose SV matrix";
			return NULL;
		}
		rhs[id] = pplhs[0];

		sr = (int)mxGetN(rhs[id]);
		sc = (int)mxGetM(rhs[id]);

		double *ptr2;
		ptr2 = mxGetPr(rhs[id]);
		ir = mxGetIr(rhs[id]);
		jc = mxGetJc(rhs[id]);

		num_samples = (int)mxGetNzmax(rhs[id]);

		elements = num_samples + sr;

		model->SV = (struct svm_node **) malloc(sr * sizeof(struct svm_node *));
		x_space = (struct svm_node *)malloc(elements * sizeof(struct svm_node));

		/* SV is in column */
		for(i=0;i<sr;i++)
		{
			int low = (int)jc[i], high = (int)jc[i+1];
			int x_index = 0;
			model->SV[i] = &x_space[low+i];
			for(j=low;j<high;j++)
			{
				model->SV[i][x_index].index = (int)ir[j] + 1; 
				model->SV[i][x_index].value = ptr2[j];
				x_index++;
			}
			model->SV[i][x_index].index = -1;
		}
#endif
		id++;
	}
#ifdef _USE_CHI_SQUARE
	if(mxIsEmpty(rhs[id]) == 0)
	{
		if(model->param.kernel_type == MULTIPLE_CHI_SQUARE)
		{
			new_featparm(model->param);
			mxGetString(rhs[id], model->param.feat_parm_text,1024);
			parse_custom_arguments(model->param.feat_parm_text, model->param.multiple_kernel_parm);
		}
	}
#endif
	mxFree(rhs);

	return model;
}
