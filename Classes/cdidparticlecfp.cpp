///////////////////////////////////////////////////////////////////////

#include "./cdidparticlecfp.h"
#include <fstream>
#include <iostream>
#include <cstdlib>

CDidParticleCFP::CDidParticleCFP()
{
  Un=-1;
}

CDidParticleCFP::~CDidParticleCFP()
{
}

void CDidParticleCFP::Clear()
{
  // delete everything
  Un=-1;
  ISF.clear();
}

void CDidParticleCFP::Insert(int nu1,int delta1,int j1,double value,
			      int delta2,int j2)
{
  // insert the isoscalara factor 'value' into the map.
  vector<int> key;
  key.push_back(nu1);
  key.push_back(delta1);
  key.push_back(j1);
  key.push_back(delta2);
  key.push_back(j2);

  // It is possible to use :
  //  ISF[key]=value;
  // but i can also give the suggested position of the new value. I assume
  // the file to be lexically ordered, thus i suggest to insert at the end.
  // This will increase speed a bit.
  pair<vector<int>,double> newpair(key,value);
  ISF.insert(ISF.end(),newpair);
}

void CDidParticleCFP::Print()
{
  vector<int> key;
  double value;
  map<vector<int>,double>::const_iterator iter;

  cout << "particle symmetry = U(" << Un << ")" << endl;
  for (iter=ISF.begin();iter!=ISF.end();iter++) {
    key=iter->first;
    value=iter->second;
    value=value*sqrt((double)key[0]); // convert isf to red.mat.
    cout << "<" << key[0] << "," << key[1] << "," << key[2] << "/2|||c+|||";
    cout << key[0]-1 << "," << key[3] << "," << key[4] << "/2> = ";
    cout << value << endl;
  }
}

double CDidParticleCFP::GetStoredISF(int nu1,int delta1,int j1,
				     int delta2,int j2)
{
  // return <nu1,delta1,j1/2|||c+|||nu1-1,delta2,j2/2>/sqrt(nu1)=ISF.
  // ATTENTION : This routine assumes that ALL ISF's which are not zero
  //             because of selection rules are stored!
  //             Every requested ISF which can not be found in the stored
  //             ones is assumed to be zero!
  double res;
  vector<int> key;
  map<vector<int>,double>::const_iterator iter;

  res=0;

  // check angular momentum selection rule
  if ((abs(j2-(Un-1))<=j1) && (j1<=j2+Un-1)) { // selection rule ok
    key.clear();
    key.push_back(nu1);
    key.push_back(delta1);
    key.push_back(j1);
    key.push_back(delta2);
    key.push_back(j2);

    iter=ISF.find(key);

    if (iter!=ISF.end()) { // isoscalar factor found
      res=iter->second;
    } else {
      printf("double CDidParticleCFP::GetStoredISF(..) :\n");
      printf("ISF not found (nu=%d). Aborting ..\n",nu1);
      exit(0);
    }
  }

  return res;
}

double CDidParticleCFP::GetCFP(int N,int nu1,int delta1,int j1,
			       int nu2,int delta2,int j2)
{
  // Return <N,nu1,delta1,j1/2||c+||N-1,nu2,delta2,j2/2>/sqrt(N).
  // It is the product of a stored isf times the isf of U(Un)>USp(Un) or
  // U(Un)>O(Un) respectively.
  // Maybe we have also to apply the reciprocal relation to the stored isf.
  // ATTENTION :If you want to get the reduced matrix element of c+, you
  // have to multiply the CFP with sqrt(N) !!!

  double res,isf1,isf2,recifactor;
  int selectionrule;
  int sil,nudiff;

  sil=Un-1;
  nudiff=nu1-nu2;

  // check some selection rules
  selectionrule=1;  // set the selection rules to be 'ok'
  if (!(N>=nu1)) { selectionrule=0; };
  if (!(N-1>=nu2)) { selectionrule=0; };
  if (!(abs(nudiff)==1)) { selectionrule=0; };
  if (!(abs(sil-j2)<=j1)) { selectionrule=0; };
  if (!(sil+j2>=j1)) { selectionrule=0; };

  //  cout << "selectionrules : " << selectionrule << endl;

  res=0;
  if (selectionrule==1) {
    // apply the reciprocal relation to 'isf1', if necessary
    if (nudiff==-1) { // we have nu1=nu2-1, apply reciprical relation
      // Get the stored isf with swithced delta1,delta2 and j1,j2
      isf1=GetStoredISF(nu2,delta2,j2,delta1,j1);
      recifactor=CalcReciprocalFactor(Un,nu2,sil,j1,j2);
      isf1=isf1*recifactor;
    } else { // no reciprocal relation
      // Get the stored isf
      isf1=GetStoredISF(nu1,delta1,j1,delta2,j2);
    }
    
    // get the isf of U(Un)>USp(Un) or U(Un)>O(Un)
    isf2=CalcISFNtoNU(Un,N,nu1,nudiff);

    //    cout << "isf1 : " << isf1 << endl;
    //    cout << "isf2 : " << isf2 << endl;

    res=isf2*isf1;
  }

  return res;
}

double CDidParticleCFP::CalcReciprocalFactor(int symmetry,int nu,int sil,
					     int j1,int j2)
{
  double res;
  int phase;
  double q1,q2;

  if (symmetry%2==0) { // fermions
    //    cout << "This is not implemented yet! Exiting .." << endl;
    //    exit(0);

    // original begin
    q1=nu*(symmetry-2*nu+4)*(j2+1);
    q2=(symmetry-nu+3)*(symmetry-2*nu+2)*(j1+1);
    //    phase=sil+j1+j2;
    // original end

    phase=j1-sil-j2;
    //    if (phase%2==0) {
    if (phase%4==0) {
      phase=-1;
    } else {
      phase=1;
    }
  } else { // bosons
    q1=nu*(2*nu+symmetry-4)*(j2+1);
    q2=(nu+symmetry-3)*(2*nu+symmetry-2)*(j1+1);
    phase=sil+j1+j2;
    //    if (phase%2==0) {
    if (phase%4==0) {
      phase=-1;
    } else {
      phase=1;
    }
  }

  //  cout << "test1"<< endl;
  //  res=CSqrtMPQ(phase,q1,q2);
  res=phase*sqrt(q1/q2);
  //  cout << "test2"<< endl;

  return res;
}

double CDidParticleCFP::CalcISFNtoNU(int symmetry,int N,int nu,int nudiff)
{
  double res;
  int phase;
  double q1,q2;

  if (symmetry%2==0) { // fermion
    if (nudiff==1) {
      phase=1;
      q1=nu*(symmetry-N-nu+2);
      q2=N*(symmetry-2*nu+2);
    } else {
      phase=-1;
      q1=(N-nu)*(symmetry-nu+2);
      q2=N*(symmetry-2*nu+2);
    }
  } else { // boson
    if (nudiff==1) {
      phase=1;
      q1=nu*(N+nu+symmetry-2);
      q2=N*(2*nu+symmetry-2);
    } else {
      phase=1;
      q1=(N-nu)*(nu+symmetry-2);
      q2=N*(2*nu+symmetry-2);
    }
  }

  //  res=CSqrtMPQ(phase,q1,q2);
  res=phase*sqrt(q1/q2);

  return res;
}

void CDidParticleCFP::MergeDFile(char* filename)
{
  // Read the isf's from a file and add them to the already stored isf's

  int nu1,delta1,j1,delta2,j2;
  double value;

  std::ifstream fin(filename, std::ios::in);
  if (!fin.is_open()) {
    std::cerr << "void CDidParticleCFP::MergeDFile(char*) :" << std::endl;
    std::cerr << "Cannot open file : " << filename << " Exiting .." << std::endl;
    std::exit(EXIT_FAILURE);
  }

  //  cout.precision(20);

  fin >> Un;
  while ((!(fin.eof())) && (!(fin.fail()))) {
    fin >> nu1;
    fin >> delta1;
    fin >> j1;
    fin >> delta2;
    fin >> j2;
    fin >> value;
    if (!(fin.fail())) {
      // cout << nu1 << " " << delta1 << " " << j1 << " " << delta2 << " ";
      // cout << j2 << " " << value << endl;
      Insert(nu1,delta1,j1,value,delta2,j2);
    }
  }

  fin.close();  
}
