//////////////////////////////////////////////////////////////////////////

#include "./csubroutines.h"
#include <chrono>

CSubroutines::CSubroutines()
{
}

CSubroutines::~CSubroutines()
{
}


// // Old implementation of SymmetricDiag, very slow for large matrices
// void CSubroutines::SymmetricDiag(vector<double>* mat,vector<double>* eigenvals,
// 		   vector<double>* eigenvecs)
// {

//     using namespace std::chrono;

//   // ---- Start timing ----
//   auto start_time = high_resolution_clock::now(); //

//   // print start message
//   cout << "Starting SymmetricDiag on matrix of size " << mat->size() << " ..." << endl;

//   // Diagonalisation of a square symmetric matrix.
//   // This routine i found in the appendix of 'The Nuclear Shell Model'
//   // of Kris Heyde

//   // We dont't need to check if the matrix size is zero, the routine
//   // handles this case.

//   if (floor(sqrt((double)mat->size()))!=ceil(sqrt((double)mat->size()))) {
//     cout << "CSubroutines::SymmetricDiag(..) :" << endl;
//     cout << "Matrix is not a square one! Exiting .." << endl;
//     exit(0);
//   }

//   int I=0,J=0,I1=0,NUL=0,NBMAX=0,IN=0,NB=0,N1=0,II=0;
//   int IDI=0,JI=0,JJ=0,K=0,IDR=0,L=0;
//   double EPS=0,EF=0,TAUSQ=0,H=0;
//   double SIG=0,COT=0,E=0,CH=0,SH=0,CO=0,SI=0,G=0,HJ=0,S1=0,S2=0;
//   double C=0,D=0,CC=0,SS=0,F=0,C1=0,C2=0,A1=0,A2=0;
//   int N;
//   int i,j,k;
//   double tmp;

//   const long MaxDim=4000;

//   N=(int)sqrt((double)mat->size());

//   if (N>MaxDim) {
//     std::cout << "SymmetricDiag: matrix dimension = " << N
//         << " (" << (N * N) << " elements, ≈ "
//         << (N * N * 8.0 / 1e6) << " MB)\n";
//     cout << "void CSubroutines::SymmetricDiag(..) :" << endl;
//     cout << "Matrix is too large. Exiting .." << endl;
//     exit(0);
//   }

//   double** A;
//   double** S;

//   A=new double*[N+1];
//   S=new double*[N+1];
//   for (i=0;i<=N;i++) {
//     A[i]=new double[N+1];
//     S[i]=new double[N+1];
//   }

//   // copy matrix to A[][]
//   for (i=0;i<N;i++) { // i is the row-index
//     for (j=0;j<N;j++) { // j is the column-index
//       A[i+1][j+1]=(*mat)[i*N+j];
//     }
//   }

//   if (N-1>=0) {
//     for (I=1;I<=N;I++) {
//       S[I][I]=1;
//       I1=I+1;
//       if (I1<=N) {
// 	for (J=I1;J<=N;J++) {
// 	  S[I][J]=0;
// 	  S[J][I]=S[I][J];
// 	}
//       }
//     }
//     NUL=0;
//     NBMAX=1000;
//     EPS=1.0E-10;
//     IN=3;
//     EF=10;
//     for (I=1;I<=IN;I++) {
//       EPS=EPS/EF;
//       for (NB=NUL;NB<=NBMAX;NB++) {
// 	if ((NB==0) || ((NB!=0) && (IDR+IDI>0))) {
// 	  IDR=0;
// 	  IDI=0;
// 	  for (II=2;II<=N;II++) {
// 	    JI=II-1;
// 	    for (JJ=1;JJ<=JI;JJ++) {
// 	      C=A[II][JJ]+A[JJ][II];
// 	      D=A[II][II]-A[JJ][JJ];
// 	      if (fabs(C)-EPS<0) {
// 		CC=1;
// 		SS=0;
// 	      } else {
// 		CC=D/C;
// 		if (CC>=0) {
// 		  SIG=1;
// 		} else {
// 		  SIG=-1;
// 		}
// 		COT=CC+SIG*sqrt(1+CC*CC);
// 		SS=SIG/sqrt(1+COT*COT);
// 		CC=SS*COT;
// 		IDR=IDR+1;
// 	      }
// 	      E=A[II][JJ]-A[JJ][II];
// 	      CH=1;
// 	      SH=0;
// 	      if (fabs(E)-EPS>=0) {
// 		CO=CC*CC-SS*SS;
// 		SI=2*SS*CC;
// 		H=0;
// 		G=0;
// 		HJ=0;
// 		for (K=1;K<=N;K++) {
// 		  if (K-II!=0) {
// 		    if (K-JJ!=0) {
// 		      H=H+A[II][K]*A[JJ][K]-
// 			A[K][II]*A[K][JJ];
// 		      S1=A[II][K]*A[II][K]+
// 			A[K][JJ]*A[K][JJ];
// 		      S2=A[JJ][K]*A[JJ][K]+
// 			A[K][II]*A[K][II];
// 		      G=G+S1+S2;
// 		      HJ=HJ+S1-S2;
// 		    }
// 		  }
// 		}
// 		D=D*CO+C*SI;
// 		H=2*H*CO-HJ*SI;
// 		F=(2*E*D-H)/(4*(E*E+D*D)+2*G);
// 		if (fabs(F)-EPS>=0) {
// 		  CH=1/sqrt(1-F*F);
// 		  SH=F*CH;
// 		  IDI=IDI+1;
// 		}
// 	      }
// 	      C1=CH*CC-SH*SS;
// 	      C2=CH*CC+SH*SS;
// 	      S1=CH*SS+SH*CC;
// 	      S2=SH*CC-CH*SS;
// 	      if (fabs(S1)+fabs(S2)!=0) {
// 		for (L=1;L<=N;L++) {
// 		  A1=A[L][II];
// 		  A2=A[L][JJ];
// 		  A[L][II]=C2*A1-S2*A2;
// 		  A[L][JJ]=C1*A2-S1*A1;
// 		  A1=S[L][II];
// 		  A2=S[L][JJ];
// 		  S[L][II]=C2*A1-S2*A2;
// 		  S[L][JJ]=C1*A2-S1*A1;
// 		}
// 		for (L=1;L<=N;L++) {
// 		  A1=A[II][L];
// 		  A2=A[JJ][L];
// 		  A[II][L]=C1*A1+S1*A2;
// 		  A[JJ][L]=C2*A2+S2*A1;
// 		}
// 	      }
// 	    }
// 	  }
// 	}
//       }
//     }
//     TAUSQ=0;
//     N1=N-1;
//     for (I=1;I<=N1;I++) {
//       II=I+1;
//       for (J=II;J<=N;J++) {
// 	TAUSQ=TAUSQ+pow((A[I][J]+A[J][I]),2);
//       }
//     }
//   } else { // N=0, do nothing
//     //    S[1][1]=1;
//   }

//   // The components of the first eigenvector are S[*][1] , the
//   // components of the i-th eigenvector are S[*][i]

//   // Sort the eigenvalues and corresponding eigenvectors, this should be
//   // replaced by a quicksort algorithm.
//   for (i=1;i<=N;i++) {
//     for (j=i+1;j<=N;j++) {
//       if (A[i][i]>A[j][j]) { // swap eigenvalue end eigenvector
// 	// swap eigenvalue
// 	tmp=A[i][i];
// 	A[i][i]=A[j][j];
// 	A[j][j]=tmp;
// 	// swap eigenvector
// 	for (k=1;k<=N;k++) {
// 	  tmp=S[k][i];
// 	  S[k][i]=S[k][j];
// 	  S[k][j]=tmp;
// 	}
//       }
//     }
//   }

//   // Copy the result to the vectors *eigenvals and *eigenvecs
//   // The first 'N' numbers in 'eigenvecs' will belong to the first eigenvector
//   eigenvals->clear();
//   eigenvecs->clear();
//   for (i=1;i<=N;i++) {
//     eigenvals->push_back(A[i][i]);
//     for (j=1;j<=N;j++) {
//       eigenvecs->push_back(S[j][i]);
//     }
//   }

//   // free the memory
//   for (i=0;i<=N;i++) {
//     delete [] A[i];
//     delete [] S[i];
//   }
//   delete [] A;
//   delete [] S;

//   // ---- End timing ----
//     auto end_time = high_resolution_clock::now();
//     auto duration_ms = duration_cast<milliseconds>(end_time - start_time).count();

//     cout << "SymmetricDiag completed in " << duration_ms << " ms (" 
//          << (duration_ms / 1000.0) << " s)" << endl;
// }

#include <iostream>
#include <vector>
#include <cmath>
#include <chrono>
#include <stdexcept>
#include <cstdlib>      // for std::getenv
#include <Eigen/Dense>

// New implementation of SymmetricDiag using Eigen (by Alex Conley)
// brew install eigen
void CSubroutines::SymmetricDiag(std::vector<double>* mat,
                                 std::vector<double>* eigenvals,
                                 std::vector<double>* eigenvecs)
{
    using namespace std;
    using namespace std::chrono;

    // Verbose controlled by env var ARB_VERBOSE=1 (or hardcode true/false)
    const bool verbose = (std::getenv("ARB_VERBOSE") && string(std::getenv("ARB_VERBOSE")) == "1");

    // Check for square matrix
    const size_t nn = mat->size();
    const int N = (int)llround(std::sqrt((double)nn));
    if ((size_t)N * (size_t)N != nn) {
        cerr << "CSubroutines::SymmetricDiag(..): matrix not square!\n";
        std::exit(1);
    }

    time_point<high_resolution_clock> t0; // declare outside the if
    if (verbose) {
        cout << "SymmetricDiag: matrix dimension = " << N
             << " (" << (size_t)N * (size_t)N << " elements, ≈ "
             << ((double)N * (double)N * 8.0 / 1e6) << " MB)\n";
        t0 = high_resolution_clock::now();
    }

    // Map the row-major flat buffer into an Eigen matrix (no copy)
    Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>
        A(mat->data(), N, N);

    // Solve: self-adjoint (symmetric) eigendecomposition
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(A);
    if (es.info() != Eigen::Success) {
        cerr << "Eigen decomposition failed\n";
        std::exit(2);
    }

    // Copy back results
    const auto& vals = es.eigenvalues();    // size N, ascending
    const auto& vecs = es.eigenvectors();   // NxN, columns = eigenvectors

    eigenvals->assign(vals.data(), vals.data() + N);
    eigenvecs->resize((size_t)N * (size_t)N);

    // Flatten eigenvectors (column j) into your row-major 1D layout
    for (int j = 0; j < N; ++j)
        for (int i = 0; i < N; ++i)
            (*eigenvecs)[(size_t)j * (size_t)N + (size_t)i] = vecs(i, j);

    if (verbose) {
        auto t1 = high_resolution_clock::now();
        auto ms = duration_cast<milliseconds>(t1 - t0).count();
        cout << "SymmetricDiag (Eigen) completed in "
             << ms << " ms (" << (ms / 1000.0) << " s)\n";
    }
}

void CSubroutines::UnfoldLocal(vector<double>* spec,vector<double>* unfspec,
			       int windowsize)
{
  // Does a kind of unfolding. The original spectrum ('spec') has
  // to hold energies which are sorted in ascending order.
  // It's the procedure of local unfolding decribed in
  // Phys. Rev. D 59 054501 in section A.2 (equations 10 and 11)
  // The parameter 'windowsize' corresponds to the parameter 'k'
  // of equation 11.

  int i,j;
  double newunfene,localmeanspacing;

  unfspec->clear();

  if (spec->size()>=1) { // take the ground state as it is
    unfspec->push_back((*spec)[0]);
  }

  for (i=1;i<(int)spec->size();i++) {
    localmeanspacing=0;
    for (j=i-windowsize;j<=i+windowsize;j++) {
      if ((j>=0) && (j+1<spec->size())) {
	localmeanspacing=localmeanspacing+(*spec)[j+1]-(*spec)[j];
      }
    }
    localmeanspacing=localmeanspacing/(2.0*windowsize+1.0);
    if (localmeanspacing==0) {
      cout << "void CSubroutines::UnfoldLocal(..) :" << endl;
      cout << "Local mean levelspacing is zero! Unable to unfold, exiting..";
      cout << endl;
      exit(0);
    }
    newunfene=(*unfspec)[i-1]+((*spec)[i]-(*spec)[i-1])/localmeanspacing;
    unfspec->push_back(newunfene);
  }
}

void CSubroutines::UnfoldGaussian(vector<double>* spec,
				  vector<double>* unfspec,double width)
{
  // Does a kind of unfolding. The original spectrum ('spec') has
  // to hold energies which are sorted in ascending order.
  // It's the procedure of Gaussian unfolding decribed in
  // Phys. Rev. D 59 054501 in section A.2 (equations 8 and 9)
  // The parameter 'width' is in units of the average level distance.

  int i,j,levelnum,intcount;
  double newunfene,groundstate,higheststate,delta,stepwidth;
  double intpos,newene;
  double opt1,chichi;

  if (width<=0) {
    cout << "void CSubroutines::UnfoldGaussian(..) :" << endl;
    cout << "The parameter 'width' is <= zero! Exiting .." << endl;
    exit(0);
  }

  unfspec->clear();

  levelnum=spec->size();
  groundstate=(*spec)[0];
  higheststate=(*spec)[spec->size()-1];

  delta=width*(higheststate-groundstate)/levelnum; // a parameter to play with
  stepwidth=delta/150.0;
  opt1=1.0/(delta*sqrt(3.1415927))*stepwidth;

  unfspec->push_back(0); // first unfolded level is zero

  newene=0;
  intpos=groundstate-2*delta;

  //  ofstream fout1("teststep.txt",ios::out);

  chichi=0;
  for (i=0;i<spec->size();i++) {
  //  for (i=1;i<100;i++) {
    // Calculate the integral (equation 8) from groundstate to (*spec)[i]
    // The integral has 'levelnum' terms (one gaussian for each level)
    //    printf("%d\n",i);
    intcount=0;
    while (intpos<=(*spec)[i]) {
      //      printf("%d\n",intcount);
      intcount=intcount+1;
      for (j=0;j<levelnum;j++) {
	newene=newene+opt1*exp(-pow((intpos-(*spec)[j])/delta,2));
      }
      //      fout1 << intpos << " " << i << " " << newene << endl;
      intpos=intpos+stepwidth;
    }
    unfspec->push_back(newene);
    chichi=chichi+pow(i-newene,2);
  }

  //  cout << chichi << endl;
  //  fout1.close();
}

vector<double> CSubroutines::CalcNNSD(vector<double>* levels,
				      double bucketsize,int bucketnum)
{
  vector<double> nnsd(bucketnum,0);
  vector<double> normspec;
  double avrgspacing,integral;
  int i,b;

  // calculate the average spacing of the spectrum 'levels[]'
  avrgspacing=0;
  for (i=0;i<levels->size()-1;i++) {
    avrgspacing=avrgspacing+((*levels)[i+1]-(*levels)[i]);
  }
  avrgspacing=avrgspacing/((double)levels->size());

  //  cout << "avrgspacing " << avrgspacing << endl;

  // get a spectrum with an average spacing of 1
  for (i=0;i<levels->size();i++) {
    normspec.push_back((*levels)[i]/avrgspacing);
  }

  // finally claculate the nnsd histogram
  for (i=0;i<normspec.size()-1;i++) {
    // calc the number of the bucket to put the spacing in
    b=(int)((normspec[i+1]-normspec[i])/bucketsize);
    // check if it is outside the arrays
    if (b<bucketnum) {  // OK, take it
      if(b<0) {
	cout << "vector<double> CSubroutines::CalcNNSD(..) :" << endl;
	cout << "Negative index! Energies are not sorted?? Exiting .." << endl;
	exit(0);
      }
      nnsd[b]=nnsd[b]+1;
    }
  }

  // normalize the distribution
  integral=0;
  for (i=0;i<bucketnum;i++) {
    integral=integral+nnsd[i];
  }
  integral=integral*bucketsize;
  if (integral==0) {
    cout << "vector<double> CSubroutines::CalcNNSD(..) :" << endl;
    cout << "NNSD is zero!" << endl;
    exit(0);
  }
  for (i=0;i<bucketnum;i++) {
    nnsd[i]=nnsd[i]/integral;
  }

  return nnsd;
}

double CSubroutines::BrodyFitNNSD(vector<double>* nnsd,double bucketsize,
				  int stepnum)
{
  int i,j;
  double chi,bestchi,omega,bestomega,tmp,x;

  if (stepnum<2) {
    cout << "double CSubroutines::BrodyFitNNSD(..) :" << endl;
    cout << "Number of steps to low! Exiting .." << endl;
    exit(0);
  }

  for (i=0;i<=stepnum-1;i++) {
    omega=(double)i/(stepnum-1);
    chi=0;
    for (j=0;j<nnsd->size();j++) {
      x=(j+0.5)*bucketsize;
      tmp=(*nnsd)[j]-(1.0+omega)*pow(x,omega)*
	exp(-pow(x,omega+1.0));
      chi=chi+tmp*tmp;
    }
    if (i==0) {
      bestchi=chi;
      bestomega=omega;
    }
    if (bestchi>chi) {
      bestchi=chi;
      bestomega=omega;
    }
  }
  return bestomega;
}

int CSubroutines::NormalizeSpectrum(vector<double>* spec,
				    vector<double>* normspec)
{
  // Take the vector 'spec' as a list of SORTED energies and
  // calculates from that a spectrum with groundstate at 0 and highest state
  // at 1. This is not always possible. If the height of the input
  // spectrum is 0 or if only one or no energies are given. In this
  // case the return value is 0, in case of success, the return value is 1.
  // If the normalization is possible, 'normspec' will receive the
  // normalized spectrum otherwise 'normspec' will be set to the
  // original spectrum.

  double shift,height;
  int i,erg;

  erg=1;
  normspec->clear();
  if (spec->size()==0) {
    erg=0;
    (*normspec)=(*spec);
  } else {
    height=(*spec)[spec->size()-1]-(*spec)[0];
    shift=(*spec)[0];
    if (height==0) {
      erg=0;
      (*normspec)=(*spec);
    } else {
      for (i=0;i<spec->size();i++) {
	normspec->push_back(((*spec)[i]-shift)/height);
      }
    }
  }

  return erg;
}

void CSubroutines::Quick_Sort(int* field,int left,int right)
{
  // quick sort of an array of integers
  // left and right specify the reagion of the array to sort (indexes)
  int pivotindex,pivotnewindex;
  if (right>left) {
    // select a pivot value
    pivotindex=left+(right-left)/2;
    pivotnewindex=QS_Partition(field,left,right,pivotindex);
    Quick_Sort(field,left,pivotnewindex-1);
    Quick_Sort(field,pivotnewindex+1,right);
  }
}

void CSubroutines::Quick_Sort_Pos(int* field,int* oriindex,int left,int right,
				  int flag)
{
  // quick sort of an array of integers and remember the original
  // positions of the elements. The original positions will be stored
  // in the (already allocated) array 'oriindex'. The first position
  // starts with zero (so its like an index).
  // left and right specify the reagion of the array to sort (indexes)
  // Set flag=0 in the inital call !!!
  int pivotindex,pivotnewindex,i;
  if (flag==0) {
    for (i=left;i<=right;i++) { oriindex[i]=i; }
  }
  // initialize the array with the original positions
  if (right>left) {
    // select a pivot value
    pivotindex=left+(right-left)/2;
    pivotnewindex=QS_Partition_Pos(field,oriindex,left,right,pivotindex);
    Quick_Sort_Pos(field,oriindex,left,pivotnewindex-1,1);
    Quick_Sort_Pos(field,oriindex,pivotnewindex+1,right,1);
  }
}

int CSubroutines::QS_Partition(int* field,int left,int right,int pivotindex)
{
  // This routine is needed by quicksort
  int pivotvalue,buf,storeindex,i;
  pivotvalue=field[pivotindex];
  // move pivot to the end (swap(field[right],field[pivotindex])
  buf=field[right];
  field[right]=field[pivotindex];
  field[pivotindex]=buf;
  storeindex=left;
  for (i=left;i<right;i++) {
    if (field[i]<=pivotvalue) {
      // swap(field[storeindex],field[i])
      buf=field[storeindex];
      field[storeindex]=field[i];
      field[i]=buf;
      storeindex=storeindex+1;
    }
  }
  // swap(field[right],field[storeindex])
  buf=field[right];
  field[right]=field[storeindex];
  field[storeindex]=buf;
  return storeindex;
}

int CSubroutines::QS_Partition_Pos(int* field,int* oriindex,
				   int left,int right,int pivotindex)
{
  // This routine is needed by quicksort
  int pivotvalue,buf,storeindex,i;
  pivotvalue=field[pivotindex];
  // move pivot to the end (swap(field[right],field[pivotindex])
  buf=field[right];
  field[right]=field[pivotindex];
  field[pivotindex]=buf;
  buf=oriindex[right];
  oriindex[right]=oriindex[pivotindex];
  oriindex[pivotindex]=buf;
  storeindex=left;
  for (i=left;i<right;i++) {
    if (field[i]<=pivotvalue) {
      // swap(field[storeindex],field[i])
      buf=field[storeindex];
      field[storeindex]=field[i];
      field[i]=buf;
      buf=oriindex[storeindex];
      oriindex[storeindex]=oriindex[i];
      oriindex[i]=buf;
      storeindex=storeindex+1;
    }
  }
  // swap(field[right],field[storeindex])
  buf=field[right];
  field[right]=field[storeindex];
  field[storeindex]=buf;
  buf=oriindex[right];
  oriindex[right]=oriindex[storeindex];
  oriindex[storeindex]=buf;
  return storeindex;
}

//////////////////////////////////////////////////////////////////////
