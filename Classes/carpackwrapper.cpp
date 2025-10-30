//////////////////////////////////////////////////////////////////////
#include "./carpackwrapper.h"
#include <fstream>   // for std::ifstream
#include <cstdlib>   // for std::exit, EXIT_FAILURE
#include <cstdio>    // for std::fprintf

CArpackWrapper::CArpackWrapper()
{
  NonZeros=0;
  NCols=0;
  NRows=0;
  IMatrix=new int[1];
  IMatrix[0]=0;
}

CArpackWrapper::~CArpackWrapper()
{
  ClearMatrix();
  delete [] IMatrix;
}

void CArpackWrapper::ClearMatrix()
{
  delete [] IMatrix;
  IMatrix=new int[1];
  IMatrix[0]=0;
  if (NonZeros>0) {
    delete [] Matrix;
    delete [] JMatrix;
  }
  NonZeros=0;
  NCols=0;
  NRows=0;
}

void CArpackWrapper::ReadMatrix(char *filename)
{
  int i;
  std::ifstream fin(filename, std::ios::in);
  if (!fin.is_open()) {
    std::fprintf(stderr, "Cannot open file: %s. Aborting..\n", filename);
    std::exit(EXIT_FAILURE);
  }

  fin >> NonZeros;
  if (NonZeros<0) {
    printf("void CArpackWrapper::ReadMatrix(char*) :\n");
    printf("Number of nonzeros<0! Aborting..\n");
    std::exit(EXIT_FAILURE);
  }

  fin >> NRows;
  if (NRows<0) {
    printf("void CArpackWrapper::ReadMatrix(char*) :\n");
    printf("Number of rows<0! Aborting..\n");
    std::exit(EXIT_FAILURE);
  }

  fin >> NCols;
  if (NCols<0) {
    printf("void CArpackWrapper::ReadMatrix(char*) :\n");
    printf("Number of columns<0! Aborting..\n");
    std::exit(EXIT_FAILURE);
  }

  if (NRows!=NCols) {
    printf("void CArpackWrapper::ReadMatrix(char*) :\n");
    printf("Matrix is not a square one! Aborting..\n");
    std::exit(EXIT_FAILURE);
  }

  IMatrix=new int[NRows+1];
  for (i=0;i<NRows+1;i++) { fin >> IMatrix[i]; };

  if (NonZeros>0) {
    JMatrix=new int[NonZeros];
    for (i=0;i<NonZeros;i++) { fin >> JMatrix[i]; };
    Matrix=new double[NonZeros];
    for (i=0;i<NonZeros;i++) { fin >> Matrix[i]; };
  }

  fin.close();
}

int CArpackWrapper::GetMatrixDim()
{
  if (NCols!=NRows) {
    printf("int CArpackWrapper::GetMatrixDim() :\n");
    printf("Matrix is not a square one. Aborting ..\n");
    exit(0);
  }
  return NRows;
}

void CArpackWrapper::Arpack_dsaupd(int eigenvalnum,double *eigenvals,
				   double *eigenvecs,int calcevecs,int whichevs,
				   int dim,int *imat,int *jmat,double *mat)
{
  // This is a wrapper routine to the arpach routine dsaupd.
  // Arguments: eigenvalnum : number of eigenvalues to calculate. The valid
  //                          range is 1 .. dim-1 (if the matrix is dim times dim).
  //            eigenvals   : pointer to an already allocated array which
  //                          will receive the eigevalues
  //            eigenvecs   : pointer to an already allocated array which
  //                          will receive the eigenvectors. Let the size of
  //                          the matrix be dim times dim, then the first dim
  //                          numbers are the components of the first eigenvector
  //                          and so on ..
  //            calcevecs   : has to be 0 or 1. For calcevecs==0 no eigenvectors
  //                          will be calculated (the argument eigenvecs does
  //                          not need to point to a valid array). For
  //                          eigenvecs==1 the eigenvectors will be calculated and
  //                          stored in the array 'eigenvecs'.
  //            whichevs    : has to be 0 or 1. For whichevs==0 the lowest
  //                          eigenvaluenum eigenvalues are calculated and for
  //                          whichevs==1 the highest ones
  //            dim,imat,jmat,mat : The matrix in ld Yale sparse matrix
  //                                format (dim by dim matrix).
  // IMPORTANT INFORMATION  : The needed storage is at maximum dim times dim
  //                          doubles , which can run you into memory problems
  //                          is the matrix is huge. The needed memory depends
  //                          on the number of requested eigenvalues !!!

  int i,j;

  // The varaibles for debug information
  int logfil,ndigit,mgetv0;
  int msaupd,msaup2,msaitr,mseigt,msapps,msgets,mseupd;
  int mnaupd,mnaup2,mnaitr,mneigh,mnapps,mngets,mneupd;
  int mcaupd,mcaup2,mcaitr,mceigh,mcapps,mcgets,mceupd;

  int maxn,maxnev,maxncv; // for the maximum dimension of the problem

  char bmat[2]="I"; /* Specifies that the right hand side matrix
		       should be the identity matrix; this makes
		       the problem a standard eigenvalue problem.
		       Setting bmat = "G" would have us solve the
		       problem Av = lBv (this would involve using
		       some other programs from BLAS, however). */

  char which[3]="__"; // this tells arpack which eigenvalues we want
  /* The possible options are
     LM: largest magnitude
     SM: smallest magnitude
     LA: largest real component
     SA: smallest real compoent
     LI: largest imaginary component
     SI: smallest imaginary component */
  // we will set this up later in the code

  if ((whichevs!=0) && (whichevs!=1)) {
    printf("void CArpackWrapper::Arpack_dsaupd(...) :\n");
    printf("Invalid parameter value given. Aborting ..\n");
    exit(0);
  }
  if ((calcevecs!=0) && (calcevecs!=1)) {
    printf("void CArpackWrapper::Arpack_dsaupd(...) :\n");
    printf("Invalid parameter value given. Aborting ..\n");
    exit(0);
  }

  maxn=2010;   // enlarge this if it is not big enougth
  //  maxnev=10;
  maxnev=maxn;
  //  maxncv=25;
  maxncv=maxn;

  if (whichevs==0) {
    which[0]='S';
    which[1]='A';
  }
  if (whichevs==1) {
    which[0]='L';
    which[1]='A';
  }

  char all[4]="All";

  int ido,n,nev,ncv,lworkl,info,ierr,ishfts,maxitr,mode1,nonv,ldv,rvec;
  double sigma;

  double tol=0.0; /* Set the tolerance. A value <= 0 specifies
		     machine recision */

  // for debug information
  ndigit=-3;
  logfil=6;
  msgets=0;
  msaitr=0;
  msapps=0;
  msaupd=1;
  msaup2=0;
  mseigt=0;
  mseupd=0;

  n=dim;  // dimension of the matrix : n by n
  ldv=n;    // leading dimension of the array v

  nev=eigenvalnum;    // number of eigenvalue requested
  ncv=4*nev;   /* The largest number of basis vectors that will
		  be used in the Implicitly Restarted Arnoldi
		  Process.  Work per major iteration is
		  proportional to N*NCV*NCV. */

  if (ncv>n) { ncv=n; }; // maximum value for ncv is n

  if (nev>=n) {
    printf("Number of requested eigenvalues has to be less than\n");
    printf("the dimension of the matrix. Aborting ..\n");
    exit(0);
  }
  if (n>maxn) {
    printf(" ERROR : n is greater than maxn \n");
    exit(0);
  }
  if (nev>maxnev) {
    printf(" ERROR : nev is greater than maxnev \n");
    exit(0);
  }
  if (ncv>maxncv) {
    printf(" ERROR : ncv is greater than maxncv \n");
    exit(0);
  }

  lworkl=ncv*(ncv+8); /* Length of the workl array */

  // now initialize arrays for input/output and working space for arpack
  double *v,*d,*workl,*workd,*resid,*ax;
  int *select,*iparam,*ipntr;
  v=new double[n*ncv];
  for (i=0;i<n*ncv;i++) { v[i]=0.0; };
  d=new double[ncv*2];
  for (i=0;i<ncv*2;i++) { d[i]=0.0; };
  workl=new double[lworkl];
  for (i=0;i<lworkl;i++) { workl[i]=0.0; };
  workd=new double[3*n];
  for (i=0;i<3*n;i++) { workd[i]=0.0; };
  resid=new double[n];
  for (i=0;i<n;i++) { resid[i]=0.0; };
  ax=new double[n];
  for (i=0;i<n;i++) { ax[i]=0.0; };
  select=new int[ncv];
  for (i=0;i<ncv;i++) { select[i]=0; };
  iparam=new int[11];
  ipntr=new int[11];
  for (i=0;i<11;i++) {
    iparam[i]=0;
    ipntr[i]=0;
  }

  info=0; /* Passes convergence information out of the iteration
	     routine. */

  // specifications for the algorithm arpack should use
  ishfts=1;
  maxitr=3*n; // maximum number of iterations
  mode1=1;
  iparam[0]=ishfts;
  iparam[2]=maxitr;
  iparam[6]=mode1;

  /* Here we enter the main loop where the calculations are
     performed.  The communication parameter ido tells us when
     the desired tolerance is reached, and at that point we exit
     and extract the solutions. */
  ido=0; // Must be zero at initial call
  do {
    dsaupd_(&ido,bmat,&n,which,&nev,&tol,resid,&ncv,v,&ldv,iparam,ipntr,
	    workd,workl,&lworkl,&info);
    if ((ido==1) || (ido==-1)) {
      /*   %--------------------------------------%
           | Perform matrix vector multiplication |
           |              y <--- OP*x             |
           | The user should supply his/her own   |
           | matrix vector multiplication routine |
           | here that takes workd(ipntr(1)) as   |
           | the input, and return the result to  |
           | workd(ipntr(2)).                     |
           %--------------------------------------% */
      // In FORTRAN the starting indext of an array is 1 but in c/c++ it is 0.
      // Therefore we have to subtract 1 from ipntr[0] and ipntr[1].
      MatrixVectorProduct(&n,imat,jmat,mat,&workd[ipntr[0]-1],&workd[ipntr[1]-1]);
    }
  } while ((ido==1) || (ido==-1));

  if (info<0) { // error
    printf("Error with _saupd, info = %d\n",info);
    printf("Check documentation in _saupd\n");
    printf("Aborting program ..\n");
    exit(0);
  } else { // np fatal error
    //    printf("No fatal error occured.\n");
    /*  %-------------------------------------------%
	| No fatal errors occurred.                 |
	| Post-Process using DSEUPD.                |
	|                                           |
	| Computed eigenvalues may be extracted.    |  
	|                                           |
	| Eigenvectors may be also computed now if  |
	| desired.  (indicated by rvec = .true.)    | 
	|                                           |
	| The routine DSEUPD now called to do this  |
	| post processing (Other modes may require  |
	| more complicated post processing than     |
	| mode1.)                                   |
	|                                           |
	%-------------------------------------------% */
    
    // extract also eigenvectors ?
    if (calcevecs==0) { rvec=0; } else { rvec=1; };
    // extract now
    dseupd_(&rvec,all,select,d,v,&ldv,&sigma,bmat,&n,which,&nev,&tol,
	    resid,&ncv,v,&ldv,iparam,ipntr,workd,workl,&lworkl,&ierr);
    if (ierr!=0) {
      printf("Error with _seupd, info = %d\n",ierr);
      printf("Check the documentation of _seupd.\n");
      printf("Aborting program ..\n");
      exit(0);
    } else {
      //      for (i=0;i<nev;i++) { printf("%f ",d[i]); };
      //      printf("\n");

      // sometimes the order of the eigenvalues is reversed, i was not able
      // to figure out when this is the case. The documentation states that
      // "the eigenvalues are given in algebraic ascending order"

      if (d[0]>d[nev-1]) { // revers order
	for (i=0;i<nev;i++) { eigenvals[i]=d[nev-1-i]; };
      } else {
	for (i=0;i<nev;i++) { eigenvals[i]=d[i]; };
      }
      if (rvec==1) { // if also eigenvecs were erquested ..
	if (d[0]>d[nev-1]) { // reverse order
	  for (i=0;i<iparam[4];i++) {
	    for (j=0;j<n;j++) {
	      eigenvecs[i*n+j]=v[(iparam[4]-1-i)*n+j];
	    }
	  }
	} else {
	  for (i=0;i<iparam[4];i++) {
	    for (j=0;j<n;j++) {
	      eigenvecs[i*n+j]=v[i*n+j];
	    }
	  }
	}
      }
    }
  }

  // clean up
  delete [] v;
  delete [] d;
  delete [] workl;
  delete [] workd;
  delete [] resid;
  delete [] ax;
  delete [] select;
  delete [] iparam;
  delete [] ipntr;
}

void CArpackWrapper::MatrixVectorProduct(int* dim,int* imat,int* jmat,double* mat,
					 double *invec,double *outvec)
{
  // Calculate the Matrix-Vector-Product with the stored Matrix in the
  // old Yale sparse matrix format.
  // Arguments : dim     : pointer to the number of columns and rows of the matrix
  //             imat,jmat and mat : the data of the matrix in old Yale
  //                                 sparse matrix format
  //             invec   : pointer to the vector which should be
  //                       multiplied by the matrix
  //             outvec  : pointer to an array which will get the result
  int i,j;
  for (i=0;i<*dim;i++) { outvec[i]=0.0; };
  for (i=0;i<*dim;i++) {
    for (j=imat[i];j<imat[i+1];j++) {
      outvec[i]=outvec[i]+mat[j]*invec[jmat[j]];
    }
  }
}
