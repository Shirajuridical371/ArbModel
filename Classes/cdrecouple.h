////////////////////////////////////////////////////////////////////////////

#if !defined _cdrecouple_h
#define _cdrecouple_h

#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <cmath>
#include <vector>
#include <map>

using namespace std;

class CDRecouple {
 public:
  vector<int> ParticleTypes; // symmetries of all different particles

  map<vector<int>,double> Termsddtt; // Map for the two body terms.
  // The key has 5 integers: n1,n2,n3,n4,n5 ; n5 in terms of hbar/2
  // [(c+_{n1}xc+_{n2})^(n5) x (c~_{n3}xc~_{n4})^(n5)]^(0)

  map<vector<int>,double> Termsdt; // Map for the one body terms.
  // The key has two integers : (c+_{n1} x c~_{n2})^(0)

  CDRecouple();
  CDRecouple(const CDRecouple&); // copy constructor
  virtual ~CDRecouple();

  CDRecouple& operator=(const CDRecouple&);

  void Add_onebody(int,int,double);
  void Add_ddtt_to_ddtt(int,int,int,int,int,double);
  void Add_dtdt_to_ddtt(int,int,int,int,int,double);

  void SetParticleTypes(vector<int>);
  void Print();
  void ClearAll();

  int GetOneBodyNum();
  int GetTwoBodyNum();
  void GetOneBodyTerms(vector<int>*,vector<double>*);
  void GetTwoBodyTerms(vector<int>*,vector<double>*);

 private:
  static const int FCT_Size=51;
  double* FCT;
  double wigner6j(int,int,int,int,int,int);
  double wigner6jhelp(int,int,int); // called by wigner6j(...)
};

#endif
