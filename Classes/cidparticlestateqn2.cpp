/////////////////////////////////////////////////////////////////

#include "./cidparticlestateqn2.h"

CidParticleStateQN2::CidParticleStateQN2()
{
  qns.clear();
  Un=-1;
  NuMax=-1;
}

CidParticleStateQN2::CidParticleStateQN2(const CidParticleStateQN2& copyqn)
{
  int i,j;
  Un=copyqn.Un;
  NuMax=copyqn.NuMax;
  if (NuMax>=0) {
    deltasizes=new int[NuMax+1];
    deltas=new int*[NuMax+1];
    for (i=0;i<=NuMax;i++) {
      deltasizes[i]=copyqn.deltasizes[i];
      deltas[i]=new int[deltasizes[i]];
      for (j=0;j<deltasizes[i];j++) { deltas[i][j]=copyqn.deltas[i][j]; };
    }
  }
  qns=copyqn.qns;
}

CidParticleStateQN2::~CidParticleStateQN2()
{
  Clear();
}

CidParticleStateQN2& CidParticleStateQN2::
operator=(const CidParticleStateQN2& assignobj)
{
  int i,j;
  if (&assignobj!=this) {
    Clear();
    Un=assignobj.Un;
    NuMax=assignobj.NuMax;
    if (NuMax>=0) {
      deltasizes=new int[NuMax+1];
      deltas=new int*[NuMax+1];
      for (i=0;i<=NuMax;i++) {
	deltasizes[i]=assignobj.deltasizes[i];
	deltas[i]=new int[deltasizes[i]];
	for (j=0;j<deltasizes[i];j++) { deltas[i][j]=assignobj.deltas[i][j]; };
      }
    }
    qns=assignobj.qns;
  }

  return (*this);
}

void CidParticleStateQN2::ConstructQN(int newn,int newun)
{
  // Arguments : 'newn'  <-> maximal particle number for which the state
  //                         quantum numbers will be constructed
  //             'newun' <-> symmetry

  ConstructQN(newn,newn,newun);
}

void CidParticleStateQN2::ConstructQN(int newn,int newmaxnu,int newun)
{
  // Arguments : 'newn'  <-> maximal particle number for which the state
  //                         quantum numbers will be constructed
  //             'newmaxnu' <-> maximum seniority of the states
  //             'newun' <-> symmetry

  int i,maxj;
  vector<int> nus,js;

  nus.clear();
  js.clear();

  // we want all spins and seniorities 0 .. newmaxnu, prepare that
  for (i=0;i<=newmaxnu;i++) { nus.push_back(i); };
  if (newun%2==0) { // fermions
    maxj=GetFermionMaxj(newun,newmaxnu);
  } else { // bosons
    maxj=newmaxnu*(newun-1);
  }
  for (i=0;i<=maxj;i++) { js.push_back(i); };

  Clear();
  CalcMultiplicities(newmaxnu,newun);
  ConstructQNfromDeltas(newn,nus,js);
}

void CidParticleStateQN2::CalcMultiplicities(int newnumax,int newun)
{
  if (newun<0) {
    printf("void CidParticleStateQN2::CalcMultiplicities(int,int) :\n");
    printf("U(n) for n<0 ?!? Aborting ..\n");
    exit(0);
  }
  if (newnumax<0) {
    printf("void CidParticleStateQN2::CalcMultiplicities(int,int) :\n");
    printf("Seniority < 0 ?!? Aborting ..\n");
    exit(0);
  }

  if (newun%2==0) { // we deal with fermions
    // for a fermion system the maximum seniority is Un/2, check this
    if (newnumax>newun/2) {
      printf("void CidParticleStateQN2::CalcMultiplicities(int,int) :\n");
      printf("Un = %d , NuMax = %d ; thats not possible, aborting ..\n",
	     newun,newnumax);
      exit(0);
    }
    CalcFermionDeltas(newnumax,newun);
  } else { // we deal with bosons
    CalcBosonDeltas(newnumax,newun);
  }
}

void CidParticleStateQN2::PrintMultiplicities()
{
  int i,j;
  for (i=0;i<=NuMax;i++) {
    for (j=0;j<deltasizes[i];j++) {
      if (Un%2==0) { // fermions
	printf("delta(%d,%d,%d/2) = %d\n",Un,i,j,deltas[i][j]);
      } else { // bosons
	printf("delta(%d,%d,%d/2) = %d\n",Un,i,j*2,deltas[i][j]);
      }
    }
  }
}

void CidParticleStateQN2::PrintAllQuantumNumbers()
{
  int n,i;
  for (n=0;n<qns.size();n++) {
    printf("N = %d :\n",n);
    for (i=0;i<qns[n].size();i++) {
      printf("|%d> = |%d,%d,%d,%d/2>\n",i,n,qns[n][i][0],
	     qns[n][i][1],qns[n][i][2]);
    }
  }
}

void CidParticleStateQN2::Clear()
{
  int i;
  // clean up
  if (NuMax>=0) {
    delete [] deltasizes;
    for (i=0;i<=NuMax;i++) { delete [] deltas[i]; };
    delete [] deltas;
  }
  Un=-1;
  NuMax=-1;
  qns.clear();
}

void CidParticleStateQN2::
ConstructQNfromDeltas(int highestn,vector<int> newsens,vector<int> newjs)
{
  int i,inu,j,ij,k,count,sil,nustart,n,jstop,ncomp;

  vector<vector<int> > emptyvecvec;
  vector<int> newqns;

  qns.clear();
  emptyvecvec.clear();

  if ((Un==-1) || (NuMax==-1)) {
    printf("void CidParticleStateQN2::ConstructQuantumNumbers(..) :\n");
    printf("You have to calculate the multiplicities first. Aborting ..");
    exit(0);
  }

  // do some checking of the parameters first
  if (Un%2==0) {
    if (highestn>Un) {
      printf("void CidParticleStateQN2::ConstructQuantumNumbers(..) :\n");
      printf("We have U(%d) particles, but Nmax = %d !? Aborting ..\n",Un,highestn);
      exit(0);
    }
  }

  for (i=0;i<newsens.size();i++) {
    if (newsens[i]>NuMax) {
      printf("void CidParticleStateQN2::ConstructQuantumNumbers(..) :\n");
      printf("Seniority out of range. Aborting ..\n");
      exit(0);
    }
  }

  sort(newsens.begin(),newsens.end());
  sort(newjs.begin(),newjs.end());

  if (Un%2==0) { // fermions
    sil=Un-1;
  } else {
    sil=(Un-1)/2;
  }

  for (n=0;n<=highestn;n++) {
    qns.push_back(emptyvecvec);
    //    printf("N = %d\n",n);
    count=1;
    for (ij=0;ij<newjs.size();ij++) { // loop for the spin
      j=newjs[ij];

      // for the given spin, calculate a lower boundary for the seniority
      if (Un%2==0) { // fermions
	nustart=0;
	i=j;
	while (i>0) {
	  i=i-((Un-1)-nustart*2);
	  nustart=nustart+1;
	}
	if (nustart%2!=n%2) { nustart=nustart+1; };
      } else { // bosons
	if (sil==0) { // s-bosons
	  nustart=n%2;
	} else {      // all other bosons
	  nustart=j/sil;
	  if (j%sil!=0) { nustart=nustart+1; };
	  if (nustart%2!=n%2) { nustart=nustart+1; };
	}
      }

      for (inu=0;inu<newsens.size();inu++) { // loop for the seniority
	i=newsens[inu];

	if (j<deltasizes[i]) {

	  // If we have fermions and a more than half filled shell,
	  // then the allowed seniority goes down.
	  ncomp=n;
	  if (Un%2==0) { if (n>Un/2) { ncomp=Un-n; }; };

	  // Check it the seniority is in the allowed range.
	  if ((i>=nustart) && (i<=ncomp) && (i%2==n%2)) {

	    for (k=0;k<deltas[i][j];k++) { // loop for the multiplicity
	      newqns.clear();
	      if (Un%2==0) { // fermions
		//		printf("|%d> = |%d,%d,%d,%d/2>\n",count,n,i,k,j);
		newqns.push_back(i);
		newqns.push_back(k);
		newqns.push_back(j);
	      } else { // bosons
		//		printf("|%d> = |%d,%d,%d,%d/2>\n",count,n,i,k,j*2);
		newqns.push_back(i);
		newqns.push_back(k);
		newqns.push_back(j*2);
	      }
	      qns[n].push_back(newqns);
	      count=count+1;
	    }

	  }
	  
	}
      }

    }
  }
}

vector<int> CidParticleStateQN2::GetStatesByQNS(int newn,vector<int> newnus,
						vector<int> newdeltas,
						vector<int> newjs)
{
  // Return a index list of states which match some quantumnumbers.
  // Arguments : newn : particle number
  //             newnu : seniorities , empty list for all
  //             newdeltas : missing labels , empty list for all
  //             newjs : spins , empty list for all

  vector<int> res;
  int i,statenum,newnusnum,newdeltasnum,newjsnum,flag,nu,delta,j,k1,k2,k3;
  int *newnusptr,*newdeltasptr,*newjsptr;

  res.clear();

  if ((newn<0) || (newn>qns.size())) {
    printf("vector<int> CidParticleQN::GetStatesByQNS(...) :\n");
    printf("I dont have states with this particle number! Exiting ..\n");
    exit(0);
  }

  newnusnum=newnus.size();
  newdeltasnum=newdeltas.size();
  newjsnum=newjs.size();
  newnusptr=&newnus.front();
  newdeltasptr=&newdeltas.front();
  newjsptr=&newjs.front();

  statenum=qns[newn].size();
  for (i=0;i<statenum;i++) {
    // if qns[n][i][0] = newnusptr[k1] for some k1 and
    // if qns[n][i][1] = newdeltasptr[k2] for some k2 and
    // if qns[n][i][2] = newjsptr[k3] for some k3 then
    // res.push_back(i)
    nu=qns[newn][i][0];
    delta=qns[newn][i][1];
    j=qns[newn][i][2];
    k1=0;
    k2=0;
    k3=0;
    flag=0;
    while ((flag==0) && (k1<newnusnum)) {
      if (nu!=newnusptr[k1]) {
	k1=k1+1;
      } else {
	flag=1;
      }
    }
    flag=0;
    while ((flag==0) && (k2<newdeltasnum)) {
      if (delta!=newdeltasptr[k2]) {
	k2=k2+1;
      } else {
	flag=1;
      }
    }
    while ((flag==0) && (k3<newjsnum)) {
      if (j!=newjsptr[k3]) {
	k3=k3+1;
      } else {
	flag=1;
      }
    }
    if (((k1<newnusnum) || (newnusnum==0)) &&
	((k2<newdeltasnum) || (newdeltasnum==0)) &&
	((k3<newjsnum) || (newjsnum==0))) {
      res.push_back(i);
    }
  }

  return res;
}

int CidParticleStateQN2::GetNu(int sn,int index)
{
  // return the spin of the (index+1)-th state with particle number 'sn'

  int res;

  if ((sn<0) || (sn>=qns.size())) {
    printf("int CidParticleQN2::GetSpin(int,int) :\n");
    printf("I dont have states with this particle number! Exiting ..\n");
    exit(0);
  }

  if ((index<0) || (index>=qns[sn].size())) {
    printf("int CidParticleQN2::GetSpin(int,int) :\n");
    printf("Index out of range! Exiting ..\n");
    exit(0);
  }

  res=qns[sn][index][0];

  return res;
}

int CidParticleStateQN2::GetDelta(int sn,int index)
{
  // return the spin of the (index+1)-th state with particle number 'sn'

  int res;

  if ((sn<0) || (sn>=qns.size())) {
    printf("int CidParticleQN2::GetSpin(int,int) :\n");
    printf("I dont have states with this particle number! Exiting ..\n");
    exit(0);
  }

  if ((index<0) || (index>=qns[sn].size())) {
    printf("int CidParticleQN2::GetSpin(int,int) :\n");
    printf("Index out of range! Exiting ..\n");
    exit(0);
  }

  res=qns[sn][index][1];

  return res;
}

int CidParticleStateQN2::GetSpin(int sn,int index)
{
  // return the spin of the (index+1)-th state with particle number 'sn'

  int res;

  if ((sn<0) || (sn>=qns.size())) {
    printf("int CidParticleQN2::GetSpin(int,int) :\n");
    printf("I dont have states with this particle number! Exiting ..\n");
    exit(0);
  }

  if ((index<0) || (index>=qns[sn].size())) {
    printf("int CidParticleQN2::GetSpin(int,int) :\n");
    printf("Index out of range! Exiting ..\n");
    exit(0);
  }

  res=qns[sn][index][2];

  return res;
}

///////////// PRIVATE FUNCTIONS /////////////

void CidParticleStateQN2::CalcFermionDeltas(int newnumax,int newun)
{
  int i,j,sil;
  int *msheme,mshemesize,maxj,**jnum;

  Clear();

  Un=newun;
  NuMax=newnumax;
  sil=Un-1; // particle angular momentum in units of hbar/2

  // allocate memory for the multiplicity and initialize it
  deltasizes=new int[NuMax+1];
  for (i=0;i<=NuMax;i++) { deltasizes[i]=GetFermionMaxj(Un,i)+1; };
  deltas=new int*[NuMax+1];
  for (i=0;i<=NuMax;i++) {
    deltas[i]=new int[deltasizes[i]];
    for (j=0;j<deltasizes[i];j++) { deltas[i][j]=0; };
  }

  // catch the NuMax=0 case
  if (NuMax==0) {
    deltas[0][0]=1;
    return;
  }

  // allocate working memory for the m-sheme and for the number of states of each j
  maxj=deltasizes[NuMax]-1;  // max ang. mom. of all states in units of hbar/2
  mshemesize=Un*(2*maxj+1); // size of the working array for the m-sheme
  msheme=new int[mshemesize];
  jnum=new int*[NuMax+1];
  for (i=0;i<=NuMax;i++) { jnum[i]=new int[deltasizes[i]]; };

  // initialize the arrays for the m-sheme
  jnum[0][0]=1; // for particle number 0 there is one state with spin 0 (vacuum)
  for (i=0;i<sil;i++) { jnum[1][i]=0; };
  jnum[1][sil]=1;
  // now pepare the working array itself
  for (i=0;i<mshemesize;i++) { msheme[i]=0; };
  // create one particle
  for (i=0;i<Un;i++) { msheme[i*(2*maxj+1)+maxj+2*i-sil]=1; };

  // all arrays are now allocate and initialized
  // do the m-sheme for seniority 2 up to NuMax amd collect after each step
  // the number of states for each spin
  for (i=2;i<=NuMax;i++) {
    AddOneFermion(Un,NuMax,msheme);
    CollectFermionJNums(Un,NuMax,i,msheme,jnum[i]);
  }

  /*
  // write some test output
  for (i=0;i<=NuMax;i++) {
    printf("N = %d\n",i);
    for (j=0;j<deltasizes[i];j++) {
      printf("J = %d/2 : %d\n",j,jnum[i][j]);
    }
  }
  */

  // The number of states for each N and l is now stored in jnum[N][l]
  // All needed information to determine the complete set of
  // quantum numbers, in particular the multiplicity is given by that array.

  // Determine now the multiplicity.
  // We store for each seniority and spin the corresponding multiplicity
  // in deltas[seniority][spin] , spin in units of hbar
  // We do this by looking at the number of states for each seniority and spin.
  // We can do seniority zero and one by hand.
  deltas[0][0]=1;
  deltas[1][sil]=1;
  // Now do seniority 2 .. NuMax
  //  printf("determining multiplicities .. ");
  for (i=2;i<=NuMax;i++) {
    for (j=0;j<deltasizes[i];j++) {
      deltas[i][j]=jnum[i][j];
      if (j<deltasizes[i-2]) { deltas[i][j]=deltas[i][j]-jnum[i-2][j]; };
    }
  }
  //  printf("done.\n");

  delete [] msheme;
  for (i=0;i<=NuMax;i++) { delete [] jnum[i]; };
  delete [] jnum;
}

int CidParticleStateQN2::GetFermionMaxj(int symm,int sen)
{
  // returns the maximum angular momentum for a system of
  // identical fermions with j=(symm-1)/2 and
  // maximum seniority 'sen' in units of hbar/2.
  // J_{max} = NuMax*j - (NuMax*NuMax-NuMax)/2 , if j is in units of hbar
  // We calculate with j in units of hbar/2, thus we do:
  // J_{max} = NuMax*j - (NuMax*NuMax-NuMax)
  int maxj;
  maxj=sen*(symm-1)-(sen*sen-sen);
  return maxj;
}

void CidParticleStateQN2::AddOneFermion(int mun,int maxsen,int *msheme)
{
  int sil,maxj,i,j,k,*nms;
  sil=mun-1;
  maxj=deltasizes[maxsen]-1;
  nms=new int[mun*(2*maxj+1)];
  for (i=0;i<mun*(2*maxj+1);i++) { nms[i]=0; };
  for (i=mun-1;i>=0;i--) {
    for (j=2*maxj+1-1;j>=0;j--) {
      if (msheme[i*(2*maxj+1)+j]!=0) {
	for (k=i+1;k<mun;k++) {
	  nms[k*(2*maxj+1)+j+2*k-sil]+=msheme[i*(2*maxj+1)+j];
	}
      }
    }
  }
  for (i=0;i<mun*(2*maxj+1);i++) { msheme[i]=nms[i]; };
  delete [] nms;
}

void CidParticleStateQN2::CollectFermionJNums(int cun,int maxsen,int sen,
					      int *msheme,int *jnum)
{
  // Collect all magnetic substates and count the number of states for each j
  // Arguments: cun    : particle symmetrie
  //            maxsen : maximum seniority (to determine the size of ms)
  //            sen    : current seniority
  //            msheme : the working array for the m-sheme
  //            jnum   : already allocated array of size=GetFermionMaxj(cun,sen)+1

  int i,maxj,maxj2,*msum;
  maxj=deltasizes[sen]-1;
  maxj2=deltasizes[maxsen]-1;
  msum=new int[maxj+1];

  // first clear the msum and lnum array
  for (i=0;i<maxj+1;i++) {
    msum[i]=0;
    jnum[i]=0;
  }

  // collect the magnetic substates
  for (i=0;i<cun*(2*maxj2+1);i++) {
    if (i%(2*maxj2+1)-maxj2>=0) {
      if (msheme[i]>0) {
	if (i%(2*maxj2+1)-maxj2>maxj) {
	  printf("void CidParticleStateQN2::CollectFermionJNums(..) :\n");
	  printf("Error during counting of states. Aborting..\n");
	  //	  printf("%d\n",i%(2*maxj2+1)-maxj2);
	  exit(0);
	}
	msum[i%(2*maxj2+1)-maxj2]+=msheme[i];
      }
    }
  }

  // now store the number of state for each j
  for (i=0;i<=maxj;i++) {
    if (i<=maxj-2) { jnum[i]=msum[i]-msum[i+2]; } else { jnum[i]=msum[i]; };
  }

  delete [] msum;
}

void CidParticleStateQN2::CalcBosonDeltas(int newnumax,int newun)
{
  int i,j,sil;
  int *msheme,mshemesize,maxl,**lnum;

  Clear();

  Un=newun;
  NuMax=newnumax;
  sil=(Un-1)/2; // particle angular momentum in units of hbar

  // allocate memory for the multiplicity and initialize it
  deltasizes=new int[NuMax+1];
  for (i=0;i<=NuMax;i++) { deltasizes[i]=i*sil+1; };
  deltas=new int*[NuMax+1];
  for (i=0;i<=NuMax;i++) {
    deltas[i]=new int[deltasizes[i]];
    for (j=0;j<deltasizes[i];j++) { deltas[i][j]=0; };
  }

  // catch the NuMax==0 case
  if (NuMax==0) {
    deltas[0][0]=1;
    return;
  }

  // allocate working memory for the m-sheme and for the number of states of each l
  maxl=deltasizes[NuMax]-1; // max ang. mom. of all states in units of hbar
  mshemesize=Un*(2*maxl+1); // size of the working array for the m-sheme
  msheme=new int[mshemesize];
  lnum=new int*[NuMax+1];
  for (i=0;i<=NuMax;i++) { lnum[i]=new int[deltasizes[i]]; };

  // initialize the arrays for the m-sheme
  lnum[0][0]=1; // for particle number 0 there is one state with spin 0 (vacuum)
  for (i=0;i<sil;i++) { lnum[1][i]=0; };
  lnum[1][sil]=1; // one particle state
  // now pepare the working array itself
  for (i=0;i<mshemesize;i++) { msheme[i]=0; };
  // create one particle
  for (i=0;i<Un;i++) { msheme[i*(2*maxl+1)+maxl+i-sil]=1; };

  // all arrays are now allocate and initialized
  // do the m-sheme for seniority 2 up to NuMax amd collect after each step
  // the number of states for each spin
  for (i=2;i<=NuMax;i++) {
    AddOneBoson(Un,NuMax,msheme);
    CollectBosonLNums(Un,NuMax,i,msheme,lnum[i]);
  }

  /*
  // write some test output
  for (i=0;i<=NuMax;i++) {
    printf("Nu = %d\n",i);
    for (j=0;j<=i*sil;j++) {
      printf("L = %d : %d\n",j,lnum[i][j]);
    }
  }
  */

  // The number of states for each N and l is now stored in lnum[N][l]
  // All needed information to determine the complete set of
  // quantum numbers, in particular the multiplicity is given by that array.

  // Determine now the multiplicity.
  // We store for each seniority and spin the corresponding multiplicity
  // in deltas[seniority][spin] , spin in units of hbar
  // We do this by looking at the number of states for each seniority and spin.
  // We can do seniority zero and one by hand.
  deltas[0][0]=1;
  deltas[1][sil]=1;
  // Now do seniority 2 .. NuMax
  //  printf("determining multiplicities .. ");
  for (i=2;i<=NuMax;i++) {
    for (j=0;j<deltasizes[i];j++) {
      deltas[i][j]=lnum[i][j];
      if (j<deltasizes[i-2]) { deltas[i][j]=deltas[i][j]-lnum[i-2][j]; };
    }
  }
  //  printf("done.\n");

  delete [] msheme;
  for (i=0;i<=NuMax;i++) { delete [] lnum[i]; };
  delete [] lnum;
}

void CidParticleStateQN2::AddOneBoson(int mun,int maxsen,int *msheme)
{
  // add one particle to the m-sheme
  int maxl,sil,i,j,k,*nms;
  sil=(mun-1)/2;
  maxl=maxsen*sil;
  nms=new int[mun*(2*maxl+1)];
  for (i=0;i<mun*(2*maxl+1);i++) { nms[i]=0; };
  for (i=mun-1;i>=0;i--) {
    for (j=2*maxl+1-1;j>=0;j--) {
      if (msheme[i*(2*maxl+1)+j]!=0) {
	for (k=i;k<mun;k++) {
	  nms[k*(2*maxl+1)+j+k-sil]+=msheme[i*(2*maxl+1)+j];
	}
      }
    }
  }
  for (i=0;i<mun*(2*maxl+1);i++) { msheme[i]=nms[i]; };
  delete [] nms;
}

void CidParticleStateQN2::CollectBosonLNums(int cun,int maxsen,int sen,
					    int *msheme,int *lnum)
{
  // Collect all magnetic substates and count the number of states for each l
  // Arguments: cun    : particle symmetrie
  //            maxsen : maximum seniority (to determine the size of ms)
  //            sen    : current seniority
  //            msheme : the working array for the m-sheme
  //            lnum   : already allocated array of size n*(Un-1)/2+1 = n*sil+1
  int i,sil,maxl,*msum;
  sil=(cun-1)/2;
  maxl=sen*sil;
  msum=new int[maxl+1];

  // first clear the msum and lnum array
  for (i=0;i<maxl+1;i++) {
    msum[i]=0;
    lnum[i]=0;
  }

  // collect the magnetic substates
  for (i=0;i<cun*(2*maxsen*sil+1);i++) {
    if (i%(2*maxsen*sil+1)-maxsen*sil>=0) {
      if (msheme[i]>0) {
	if (i%(2*maxsen*sil+1)-maxsen*sil>maxl) {
	  printf("void CidParticleStateQN2::CollectBosonLNums(..) :\n");
	  printf("Error during counting of states. Aborting..\n");
	  //	  printf("%d\n",i%(2*maxsen*sil+1)-maxsen*sil);
	  exit(0);
	}
	msum[i%(2*maxsen*sil+1)-maxsen*sil]+=msheme[i];
      }
    }
  }

  // now store the number of state for each l
  for (i=0;i<=maxl;i++) {
    if (i!=maxl) { lnum[i]=msum[i]-msum[i+1]; } else { lnum[i]=1; };
  }

  delete [] msum;
}
