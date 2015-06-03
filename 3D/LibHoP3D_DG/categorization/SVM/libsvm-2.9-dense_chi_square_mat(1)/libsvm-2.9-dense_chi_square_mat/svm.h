#ifndef _LIBSVM_H
#define _LIBSVM_H

#define LIBSVM_VERSION 290

#ifdef _ALL_FLOAT
typedef float SVM_REAL;
#else
typedef double SVM_REAL;
#endif

#ifdef __cplusplus
extern "C" {
#endif

extern int libsvm_version;

#ifdef _DENSE_REP
struct svm_node
{
	int dim;
	SVM_REAL *values;
};

struct svm_problem
{
	int l;
	SVM_REAL *y;
	struct svm_node *x;
};

#else
struct svm_node
{
	int index;
	SVM_REAL value;
};

struct svm_problem
{
	int l;
	SVM_REAL *y;
	struct svm_node **x;
};
#endif

enum { C_SVC, NU_SVC, ONE_CLASS, EPSILON_SVR, NU_SVR };	/* svm_type */
#ifdef _USE_CHI_SQUARE
enum { LINEAR, POLY, RBF, SIGMOID, PRECOMPUTED, CHI_SQUARE, MULTIPLE_CHI_SQUARE }; /* kernel_type */
#else
enum { LINEAR, POLY, RBF, SIGMOID, PRECOMPUTED }; /* kernel_type */
#endif

#ifdef _USE_CHI_SQUARE
struct FEAT_PARM
{
    int length;
    float width;
    float weight;
    struct FEAT_PARM *next;
};
#endif

struct svm_parameter
{
	int svm_type;
	int kernel_type;
	int degree;	/* for poly */
	SVM_REAL gamma;	/* for poly/rbf/sigmoid */
	SVM_REAL coef0;	/* for poly/sigmoid */

	/* these are for training only */
	SVM_REAL cache_size; /* in MB */
	SVM_REAL eps;	/* stopping criteria */
	SVM_REAL C;	/* for C_SVC, EPSILON_SVR and NU_SVR */
	int nr_weight;		/* for C_SVC */
	int *weight_label;	/* for C_SVC */
	SVM_REAL* weight;		/* for C_SVC */
	SVM_REAL nu;	/* for NU_SVC, ONE_CLASS, and NU_SVR */
	SVM_REAL p;	/* for EPSILON_SVR */
	int shrinking;	/* use the shrinking heuristics */
	int probability; /* do probability estimates */
#ifdef _USE_CHI_SQUARE
	FEAT_PARM *multiple_kernel_parm;
	char feat_parm_text[1024];
#endif
};

#ifdef _USE_CHI_SQUARE

void free_featparm(FEAT_PARM *feats);
void parse_custom_arguments(char *argument, FEAT_PARM *feats);
SVM_REAL kernel_multiple_chi_square(int min_dim, const SVM_REAL *x, const SVM_REAL *y, const FEAT_PARM * multiple_kernel_parm);
void new_featparm(svm_parameter &param);
#endif

//
// svm_model
//
struct svm_model
{
	svm_parameter param;	// parameter
	int nr_class;		// number of classes, = 2 in regression/one class svm
	int l;			// total #SV
#ifdef _DENSE_REP
	svm_node *SV;		// SVs (SV[l])
#else
	svm_node **SV;		// SVs (SV[l])
#endif
	SVM_REAL **sv_coef;	// coefficients for SVs in decision functions (sv_coef[k-1][l])
	SVM_REAL *rho;		// constants in decision functions (rho[k*(k-1)/2])
	SVM_REAL *probA;		// pariwise probability information
	SVM_REAL *probB;

	// for classification only

	int *label;		// label of each class (label[k])
	int *nSV;		// number of SVs for each class (nSV[k])
				// nSV[0] + nSV[1] + ... + nSV[k-1] = l
	// XXX
	int free_sv;		// 1 if svm_model is created by svm_load_model
				// 0 if svm_model is created by svm_train
};

struct svm_model *svm_train(const struct svm_problem *prob, const struct svm_parameter *param);
void svm_cross_validation(const struct svm_problem *prob, const struct svm_parameter *param, int nr_fold, SVM_REAL *target);

int svm_save_model(const char *model_file_name, const struct svm_model *model);
struct svm_model *svm_load_model(const char *model_file_name);

int svm_get_svm_type(const struct svm_model *model);
int svm_get_nr_class(const struct svm_model *model);
void svm_get_labels(const struct svm_model *model, int *label);
SVM_REAL svm_get_svr_probability(const struct svm_model *model);

void svm_predict_values(const struct svm_model *model, const struct svm_node *x, SVM_REAL* dec_values);
SVM_REAL svm_predict(const struct svm_model *model, const struct svm_node *x);
SVM_REAL svm_predict_probability(const struct svm_model *model, const struct svm_node *x, SVM_REAL* prob_estimates);

void svm_destroy_model(struct svm_model *model);
void svm_destroy_param(struct svm_parameter *param);

const char *svm_check_parameter(const struct svm_problem *prob, const struct svm_parameter *param);
int svm_check_probability_model(const struct svm_model *model);

extern void (*svm_print_string) (const char *);

#ifdef __cplusplus
}
#endif

#endif /* _LIBSVM_H */
