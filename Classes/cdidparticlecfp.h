///////////////////////////////////////////////////////////////////////

#if !defined _cdidparticlecfp_h
#define _cdidparticlecfp_h

#include <fstream>
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <map>

using namespace std;

class CDidParticleCFP {
  // Class for reading ISF O(Un) > O(3) or USp(Un) > SU(2) resp. from a file
  // and storing them in double precision. This class can then serve as a
  // database in some model code.
 public:
  int Un; // symmetry of the particle, Un=5 for d-boson, Un=4 for 3/2-fermions

  // This map contains <nu,delta1,j1|||c+|||nu-1,delta2,j2>/sqrt(nu)=ISF
  map<vector<int>,double> ISF; // key - values pairs for the ISF's

  CDidParticleCFP();
  virtual ~CDidParticleCFP();

  void Clear();
  void Insert(int,int,int,double,int,int);
  void Print();

  double GetStoredISF(int,int,int,int,int); // get a stored isf
  double GetCFP(int,int,int,int,int,int,int); // get a CFP

  double CalcReciprocalFactor(int,int,int,int,int);
  double CalcISFNtoNU(int,int,int,int);

  void MergeDFile(char* filename); // merge the ISF's from file into database
};

#endif
