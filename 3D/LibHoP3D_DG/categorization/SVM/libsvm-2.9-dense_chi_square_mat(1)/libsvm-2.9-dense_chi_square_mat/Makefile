MATLABDIR ?= /usr/local/MATLAB/R2010b

CXX ?= g++
CC = gcc
CFLAGS = -Wall -Wconversion -O3 -fPIC -march=nocona -D _DENSE_REP -D _USE_CHI_SQUARE -I$(MATLABDIR)/extern/include

OMPFLAGS = -fopenmp

FLOATFLAGS = -D _ALL_FLOAT

MEX = $(MATLABDIR)/bin/mex
MEX_OPTION = CC\#$(CXX) CXX\#$(CXX) CFLAGS\#"$(CFLAGS) $(OMPFLAGS)" CXXFLAGS\#"$(CFLAGS) $(OMPFLAGS)"
MEX_OPTION += -largeArrayDims
MEX_EXT = $(shell $(MATLABDIR)/bin/mexext)
LDIRS= -L$(MATLABDIR)/bin/glnxa64
MEX_LIBS = -lmex -shared

CFLAGS_CHI = -Wall -Wconversion -O3 -fPIC -march=nocona -ffast-math -fomit-frame-pointer
SHVER = 1

all: binary binary_float

float: binary_float

binary: svmpredict_chi2.$(MEX_EXT) svmtrain_chi2.$(MEX_EXT) libsvmread.$(MEX_EXT) libsvmwrite.$(MEX_EXT)

binary_float: svmpredict_chi2_float.$(MEX_EXT) svmtrain_chi2_float.$(MEX_EXT) libsvmread.$(MEX_EXT) libsvmwrite.$(MEX_EXT)

svmpredict_chi2_float.$(MEX_EXT):     svmpredict.c svm.h svm_float.o svm_model_matlab_float.o chi2float.o
	$(CXX) $(LDIRS) $(CFLAGS) $(FLOATFLAGS) $(OMPFLAGS) $(MEX_LIBS) svmpredict.c svm_float.o svm_model_matlab_float.o chi2float.o -o svmpredict_chi2_float.$(MEX_EXT)

svmtrain_chi2_float.$(MEX_EXT):       svmtrain.c svm.h svm_float.o svm_model_matlab_float.o chi2float.o
	$(CXX) $(LDIRS) $(CFLAGS) $(FLOATFLAGS) $(OMPFLAGS) $(MEX_LIBS) svmtrain.c svm_float.o svm_model_matlab_float.o chi2float.o -o svmtrain_chi2_float.$(MEX_EXT)

svmpredict_chi2.$(MEX_EXT):     svmpredict.c svm.h svm.o svm_model_matlab.o chi2double.o
	$(CXX) $(LDIRS) $(CFLAGS) $(OMPFLAGS) $(MEX_LIBS) svmpredict.c svm.o svm_model_matlab.o chi2double.o -o svmpredict_chi2.$(MEX_EXT)

svmtrain_chi2.$(MEX_EXT):       svmtrain.c svm.h svm.o svm_model_matlab.o chi2double.o
	$(CXX) $(LDIRS) $(CFLAGS) $(OMPFLAGS) $(MEX_LIBS) svmtrain.c svm.o svm_model_matlab.o chi2double.o -o svmtrain_chi2.$(MEX_EXT)

libsvmread.$(MEX_EXT):	libsvmread.c
	$(MEX) $(MEX_OPTION) libsvmread.c

libsvmwrite.$(MEX_EXT):	libsvmwrite.c
	$(MEX) $(MEX_OPTION) libsvmwrite.c

svm_model_matlab.o:     svm_model_matlab.c svm.h
	$(CXX) $(CFLAGS) -c svm_model_matlab.c
	
svm_model_matlab_float.o: svm_model_matlab.c svm.h
	$(CXX) $(CFLAGS) $(FLOATFLAGS) -o svm_model_matlab_float.o -c svm_model_matlab.c

svm.o: svm.cpp svm.h
	$(CXX) $(CFLAGS) -c $(OMPFLAGS) -o svm.o svm.cpp
svm_float.o: svm.cpp svm.h
	$(CXX) $(CFLAGS) $(FLOATFLAGS) -c $(OMPFLAGS) -o svm_float.o svm.cpp
chi2double.o : chi2double.c chi2double.h Makefile
	$(CXX) $(CFLAGS_CHI) -c $(OMPFLAGS) -o chi2double.o chi2double.c
chi2float.o : chi2float.c chi2float.h Makefile
	$(CXX) $(CFLAGS_CHI) -c $(OMPFLAGS) -o chi2float.o chi2float.c

clean:
	rm -f *~ *.o *.$(MEX_EXT) *.obj *.mex*
