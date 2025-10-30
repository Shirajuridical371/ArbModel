/////////////////////////////////////////////////////////////////////////////

#if !defined _carpackwrapper_h
#define _carpackwrapper_h

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>

using namespace std;

// Decaration of ARPACK routines
extern "C" {
  void dsaupd_(const int *ido,const char *bmat,const int *n,
	       const char *which,const int *nev,const double *tol,
	       const double *resid,const int *ncv,const double *v,
	       const int *ldv,const int *iparam,const int *ipntr,
	       const double *workd,const double *workl,const int *lworkl,
	       const int *info);

  void dseupd_(const int *rvec,const char *all,const int *select,
	       const double *d,const double *z,const int *ldz,
	       const double *sigma,const char *bmat,const int *n,
	       const char *which,const int *nev,const double *tol,
	       const double *resid,const int *ncv,const double *v,
	       const int *ldv,const int *iparam,const int *ipntr,
	       const double *workd,const double *workl,const int *lworkl,
	       const int *ierr);
}

class CArpackWrapper {
 public:
  CArpackWrapper();
  virtual ~CArpackWrapper();

  void ClearMatrix();
  void ReadMatrix(char*);
  int GetMatrixDim();

  void Arpack_dsaupd(int,double*,double*,int,int,int,int*,int*,double*);

  void MatrixVectorProduct(int*,int*,int*,double*,double*,double*);

  // for the storage of a matrix in the old yale sparse matrix format
  int NRows,NCols,NonZeros;
  int *IMatrix,*JMatrix;
  double *Matrix;
};

#endif
