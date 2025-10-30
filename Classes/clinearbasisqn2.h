///////////////////////////////////////////////////////////////////////

#if !defined _clinearbasis2_h
#define _clinearbasis2_h

#include <cstdlib>
#include <iostream>
#include <vector>

#include "./cidparticlestateqn2.h"

using namespace std;

class CLinearBasisQN2 {
 public:
  CLinearBasisQN2();
  CLinearBasisQN2(const CLinearBasisQN2&); // copy constructor
  virtual ~CLinearBasisQN2();

  vector<int> Basis; // The basis in terms of quantum numbers
  // Format : Each state is represented by ptypes.size()*5 quantum numbers.
  // For each particle type we have : N,nu,delta,j and the again an angular
  // momentum which is the coupled spin of the identical particle state
  // (which is represented by N,nu,delta,j) and the intermediate spin
  // of the coupling before (which is zero if we deal with the first one).

  vector<int> ParticleConfigs;    // These two arrays contain all configurations
  vector<int> ConfigStateIndexes; // of particles and for each configuration the
  // index of the first and last state as well as the index of the first and last
  // combination which belongs to that configuration.
  // Each configuration consists of ptypes.size() numbers.

  vector<int> IdPartStateCombis; // These two arrays contain all combination
  vector<int> CombiStateIndexes; // of indexes of identical particle states and
  // the indexes of the first and last state correspondig to that combination.
  // The particle numbers are not included in the combination, just state indexes
  // corresponding to the proper particle numbers.

  vector<int> ptypes;   // symmetries of the different particles
  vector<int> maxns;    // max. particle number of the different particles
  vector<int> limits;   // restrictions of the basis
  vector<int> parities; // parities of the different particles
  int N;                // all over all particle number of the states

  CLinearBasisQN2& operator=(const CLinearBasisQN2&);

  void ConstructBasis(int,int,int,int,vector<int>,vector<int>,
		      vector<int>,vector<int>);
  void ConstructBasis(int,int,int,int,vector<int>,vector<int>,
		      vector<int>,vector<int>,vector<int>);
  void ConstructBasis(int,int,int,int,vector<int>,vector<int>,
		      vector<int>,vector<int>,vector<vector<int> >);
  void GetAllCouplingShemes(int*,int,int);
  int GetNextSheme(int*,int*,int,int);
  void StoreStateQN(int*,int*,int);

  //  void FilterSeniority(vector<int>); // filter out states with some seniority

  void Print();
  int GetStateNum();
  int GetConfNum();
  void Clear();
  void PrintConfigs();
  void PrintCombis();
  int* GetPtr(int);

  void GetJArray(vector<int>*);

  void SelectStates(int*,int,vector<int>*);
  int GetConfIndex(int*);

 private:
  void PrintSheme(int*,int*,int);
  int FindConfigurations(int,int,int*,int*,int*,int,
			 int,int*,int,int*,int*);
  int LimitsOK(int*,int*,int,int,int);
  int GetParity(int*,int*,int);

  ////////////////////////////////////////////////////////////////////
  // The following global data is needed for the construction process
  int* configurations;
  int confignum,configindex;

  vector<CidParticleStateQN2> idpartqn;

  vector<vector<int> > indexlists;
  int* indexes;
  int J;
  //////////////////////////////////////////////
};

#endif
