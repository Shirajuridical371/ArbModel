///////////////////////////////////////////////////////////////////////

#if !defined _cdreducedmatrix2_h
#define _cdreducedmatrix2_h

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <iostream>
#include <fstream>

#include "./clinearbasisqn2.h"
#include "./cdidparticlecfp.h"

using namespace std;

class CDReducedMatrix2 {
  // Class for constructing and using sparse matrices consisting of
  // reduced matrix elements of operators (Talmi convention).
  // <a1,j1,m1|O^j_m|a2,j2,m2> = (-)^(j1-m1) * wigner3j(j1,j,j2,-m1,m,m2) *
  //                             <a1,j1||O^j||a2,j2>

  // The sparse matrix is stored in the OLD YALE SPARSE MATRIX FORMAT:
  // Let the matrix have 'NRows' rows and 'NCols' columns.
  // The nonzero elements of the matrix are stored in an array 'Matrix' in row
  // format. Let the number of nonzero entries be 'NonZeros'. Then 'Matrix'
  // has the 'NonZeros' elements.
  // There are two more arrays 'IMatrix' and 'JMatrix'.
  // 'IMatrix' has 'NRows'+1 elements. The first
  // element of row i is 'Matrix[IMatrix[i]]'. The number of nonzero elements
  // of row i is 'IMatrix[i+1]'-'Matrix[i]', which is why 'IMatrix' has 'NRows+1'
  // elements. The column index of the element 'Matrix[i]' is 'JMatrix[i]'.
  // So the size of 'JMatrix' is 'NonZeros'.
  // IMPORTANT NOTE: The column numbers in JMatrix for a given row a NOT
  // always in ascending order!

 public:
  CDReducedMatrix2();
  CDReducedMatrix2(const CDReducedMatrix2&);  // copy constructor
  virtual ~CDReducedMatrix2();

  CDReducedMatrix2& operator=(const CDReducedMatrix2&);
  CDReducedMatrix2 operator+(CDReducedMatrix2);
  CDReducedMatrix2 operator*(double);

  void Mul(double);
  void MulVec(double*,double*); // matrix vector multiplication
  void Mul3jAndPhase(int,int,int,int,int,int,int);
  void Transpose(CDReducedMatrix2*);

  double GetRMatElement(int*,int*,int,int,int*,CDidParticleCFP*);

  void CalcReducedMatrix(CLinearBasisQN2*,CLinearBasisQN2*,
			 CDidParticleCFP*,int);
  void CalcReducedMatrix(CDReducedMatrix2*,CDReducedMatrix2*,
			 int*,int*,int*,int,int,int);

  void ClearMatrix();
  void SetNullMatrix(int,int);
  void Print();
  int GetRowNum();
  int GetColumnNum();
  void GetMatrix(double*);
  void WriteToFile(char*);

  // The matrix itself:
  int NonZeros,NRows,NCols;
  double *Matrix;
  int *IMatrix,*JMatrix;

 private:
  // for the calculation of wigner3j- and wigner6j-symbols:
  double wigner6jhelp(int,int,int);
  double wigner6j(int,int,int,int,int,int);
  double wigner3j(int,int,int,int,int,int);
  double clebsch(int,int,int,int,int,int);
  double* FCT;
  static const int FCT_Size=71;
};

#endif
