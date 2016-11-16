#include "mex.h"
#include "math.h"
#define EPS  2.2204460492503131e-016
#include <string.h>
#include "multiple_os_thread.h"
double pow2(double a){return a*a;} 

voidthread chifunction(double **Args) {
    int i, j, k;
    int indexC, offset1, offset2;
    double *Nthreadsd;
    int Nthreads;
    double *F1, *F2;
    double *P1, *P2;
    double *C;
    double *MD;
    int F1size[2], F2size[2];
    int P1size[2], P2size[2];
    int Csize[2];
    double *F1sized, *F2sized;
    double *P1sized, *P2sized;
    double *Csized;
    double f1, f2;
    /* offset */
    int ThreadOffset=0;
    double maxdist2;
    /* The thread ID number/name */
    double *ThreadID;
    double dist;
    /* Get the inputs */
    F1=Args[0];
    F1sized=Args[1];
    F2=Args[2];
    F2sized=Args[3];
    P1=Args[4];
    P1sized=Args[5];
    P2=Args[6];
    P2sized=Args[7];
    MD=Args[8];
    C=Args[9];
    Csized=Args[10];
    ThreadID=Args[11];
    Nthreadsd=Args[12];
        
    /* Maximum quadratic matching distance */
    maxdist2=MD[0]*MD[0];
    
    /* Sizes to int */
    F1size[0]=(int)F1sized[0]; F1size[1]=(int)F1sized[1];
    F2size[0]=(int)F2sized[0]; F2size[1]=(int)F2sized[1];
    P1size[0]=(int)P1sized[0]; P1size[1]=(int)P1sized[1];
    P2size[0]=(int)P2sized[0]; P2size[1]=(int)P2sized[1];
    Csize[0] =(int)Csized[0];  Csize[1] =(int)Csized[1];
    

    Nthreads=(int)Nthreadsd[0];
    ThreadOffset=(int) ThreadID[0];
    
    /* Calculate the CHI-square distance between feature-vectors*/
    for(i=ThreadOffset; i<F2size[1]; i+=Nthreads) {
        offset2=i*F2size[0];
        for(k=0; k<F1size[1]; k++) {
            indexC=k+i*Csize[0];
            offset1=k*F1size[0];
            if(P1size[1]>2)
            {
                dist=pow2(P2[i]-P1[k])+pow2(P2[i+P2size[0]]-P1[k+P1size[0]])+pow2(P2[i+2*P2size[0]]-P1[k+2*P1size[0]]);
            }
            else
            {
                dist=pow2(P2[i]-P1[k])+pow2(P2[i+P2size[0]]-P1[k+P1size[0]]);
            }    
            
            if(dist>maxdist2)
            {
                C[indexC]=1e100;
            }
            else
            {
                for(j=0; j<F1size[0]; j++) {
                    f1=F1[offset1+j];
                    f2=F2[offset2+j];
                    C[indexC]+=((f1-f2)*(f1-f2))/(f1+f2+EPS);
                }
            }
        }
    }
    /*  explicit end thread, helps to ensure proper recovery of resources allocated for the thread */
    EndThread;
}

/* The matlab mex function */
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {
    double *F1, *F2;
    double *P1, *P2;
    double *C;
    double *MD;
    const mwSize *F1dims, *F2dims;
    const mwSize *P1dims, *P2dims;
    int F1size[2], F2size[2];
    int P1size[2], P2size[2];
    int Csize[2];
    double F1sized[2], F2sized[2];
    double P1sized[2], P2sized[2];
    double Csized[2];
    int Nthreads;
    double Nthreadsd[1];
	int i;
    /* double pointer array to store all needed function variables) */
    double ***ThreadArgs;
    double **ThreadArgs1;
    /* ID of Threads */
    double **ThreadID;              
    double *ThreadID1;
	/* Handles to the worker threads */
    ThreadHANDLE *ThreadList;
        
    /* Get the input feature-matrices */
    F1=mxGetPr(prhs[0]);
    F2=mxGetPr(prhs[1]);
    P1=mxGetPr(prhs[2]);
    P2=mxGetPr(prhs[3]);
    MD=mxGetPr(prhs[4]);
    
    /* Get input dimensions*/
    F1dims = mxGetDimensions(prhs[0]);
    F2dims = mxGetDimensions(prhs[1]);
    P1dims = mxGetDimensions(prhs[2]);
    P2dims = mxGetDimensions(prhs[3]);
    
    /* Convert input dimensions from const mwsize to int */
    F1size[0]=F1dims[0]; F1size[1]=F1dims[1];
    F2size[0]=F2dims[0]; F2size[1]=F2dims[1];
    P1size[0]=P1dims[0]; P1size[1]=P1dims[1];
    P2size[0]=P2dims[0]; P2size[1]=P2dims[1];
    
    Csize[0]=F1size[1]; Csize[1]=F2size[1];

    /* Create the ouput Cost matrix */
    plhs[0] = mxCreateNumericArray(2, Csize, mxDOUBLE_CLASS, mxREAL);
    C=mxGetPr(plhs[0]);
    memset(C, 0, sizeof(C));
    
    /* Sizes to double, for muli-threading*/
    F1sized[0]=F1size[0]; F1sized[1]=F1size[1];
    F2sized[0]=F2size[0]; F2sized[1]=F2size[1];
    P1sized[0]=P1size[0]; P1sized[1]=P1size[1];
    P2sized[0]=P2size[0]; P2sized[1]=P2size[1];
    Csized[0] =Csize[0];  Csized[1] = Csize[1];
    
    /* Get max number of threads */
    Nthreads=getNumCores();
    Nthreadsd[0]=Nthreads;

    /* Reserve room for handles of threads in ThreadList  */
    ThreadList = (ThreadHANDLE*)malloc(Nthreads* sizeof( ThreadHANDLE ));
    ThreadID = (double **)malloc( Nthreads* sizeof(double *) );
    ThreadArgs = (double ***)malloc( Nthreads* sizeof(double **) );
    
    for (i=0; i<Nthreads; i++) {
        /*  Make Thread ID  */
        ThreadID1= (double *)malloc( 1* sizeof(double) );
        ThreadID1[0]=i;
        ThreadID[i]=ThreadID1;
        
        /*  Make Thread Structure  */
        ThreadArgs1 = (double **)malloc( 13* sizeof( double * ) );
        ThreadArgs1[0]=F1;
        ThreadArgs1[1]=F1sized;
        ThreadArgs1[2]=F2;
        ThreadArgs1[3]=F2sized;
        ThreadArgs1[4]=P1;
        ThreadArgs1[5]=P1sized;
        ThreadArgs1[6]=P2;
        ThreadArgs1[7]=P2sized;
        ThreadArgs1[8]=MD;
        ThreadArgs1[9]=C;
        ThreadArgs1[10]=Csized;
        ThreadArgs1[11]=ThreadID[i];
        ThreadArgs1[12]=Nthreadsd;
        
        /* Start a Thread  */
        ThreadArgs[i]=ThreadArgs1;
        StartThread(ThreadList[i], &chifunction, ThreadArgs[i])
    }
    
    /* Wait for all threads to finish */
    for (i=0; i<Nthreads; i++) { WaitForThreadFinish(ThreadList[i]); }
    
    /* Remove all threading resources */
    for (i=0; i<Nthreads; i++) { free(ThreadArgs[i]); free(ThreadID[i]); }
    free(ThreadArgs); free(ThreadID); free(ThreadList);
}

