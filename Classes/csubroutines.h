//////////////////////////////////////////////////////////////////////

#if !defined _csubroutines_h
#define _csubroutines_h

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <cmath>
#include <vector>

using namespace std;

class CSubroutines {
 public:
  CSubroutines();
  virtual ~CSubroutines();

  void SymmetricDiag(vector<double>*,vector<double>*,vector<double>*);

  void UnfoldLocal(vector<double>*,vector<double>*,int);
  void UnfoldGaussian(vector<double>*,vector<double>*,double);

  vector<double> CalcNNSD(vector<double>*,double,int);

  // args: nnsd,bucketsize,steps from 0 and 1 for omega (101 for 0.0, 0.1, ..)
  double BrodyFitNNSD(vector<double>*,double,int);

  int NormalizeSpectrum(vector<double>*,vector<double>*);

  void Quick_Sort(int*,int,int);  // quick sort integers

  // quick sort+remember ori. pos
  void Quick_Sort_Pos(int*,int*,int,int,int flag=0); 

 private:
  int QS_Partition(int*,int,int,int); // needed by Quick_Sort
  int QS_Partition_Pos(int*,int*,int,int,int); // needed by Quick_Sort_Pos
};

#endif
