///////////////////////////////////////////////////////////////////////

#if !defined _cidparticlestateqn2_h
#define _cidparticlestateqn2_h

#include <cstdlib>
#include <cstdio>

#include <vector>
#include <algorithm>

using namespace std;

class CidParticleStateQN2 {
  // This class can be used for two things:
  // 1. Getting quantum numbers of states of identical particles
  //    (bosons and fermions of arbitrary spin). The group chain is:
  //    for bosons   : U(n) -  SO(n) - SO(3)
  //    for fermions : U(n) - USp(n) - SU(2)
  // 2. Calculating the multiplicity between SO(n) - SO(3) (for bosons) or
  //    USp(n) - SU(2) (for fermions)
 public:
  // qns[n][i][0] is seniority if the (i+1)-th state with N = n
  // qns[n][i][1] is missing label if the (i+1)-th state with N = n
  // qns[n][i][2] is spin in units of hbar/2 if the (i+1)-th state with N = n
  vector<vector<vector<int> > > qns; // all state quantum numbers

  CidParticleStateQN2();  // standard constructor
  CidParticleStateQN2(const CidParticleStateQN2&); // copy constructor
  virtual ~CidParticleStateQN2();

  CidParticleStateQN2& operator=(const CidParticleStateQN2&);

  void ConstructQN(int,int); // for N=0 up to something
  void ConstructQN(int,int,int); // with some given maximum seniority

  void CalcMultiplicities(int,int);
  void PrintMultiplicities();
  void PrintAllQuantumNumbers();
  void Clear();
  void ConstructQNfromDeltas(int,vector<int>,vector<int>);

  vector<int> GetStatesByQNS(int,vector<int>,vector<int>,vector<int>);

  int GetNu(int,int);
  int GetDelta(int,int);
  int GetSpin(int,int);

 private:
  // For the fermion multiplicities
  void CalcFermionDeltas(int,int);
  int GetFermionMaxj(int,int);
  void AddOneFermion(int,int,int*);
  void CollectFermionJNums(int,int,int,int*,int*);

  // For the boson multiplicities
  void CalcBosonDeltas(int,int);
  void AddOneBoson(int,int,int*);
  void CollectBosonLNums(int,int,int,int*,int*);

  // The Array of arrays 'deltas' holds the multiplicities.
  // deltas[nu][j] is the multiplicity corresponding to seniority nu and spin j
  // If Un%2==0 then j is in units of hbar/2, if Un%2==1 then j is in units of hbar.
  int **deltas;    // (NuMax+1) arrays of sizes deltasizes[i] with i=0 .. NuMax
  int *deltasizes; // for storing the (NuMax+1) sizes of deltas[i] with i=0 .. NuMax
  int Un;          // symmetry of the particles , -1 for 'nothing is stored'
  int NuMax;       // maximum seniority of the deltas, -1 for 'nothing is stored'
};

#endif
