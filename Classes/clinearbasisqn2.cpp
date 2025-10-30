///////////////////////////////////////////////////////////////////////

#include "./clinearbasisqn2.h"

CLinearBasisQN2::CLinearBasisQN2()
{
  Basis.clear();
  ptypes.clear();
  maxns.clear();
  limits.clear();
  ParticleConfigs.clear();
  ConfigStateIndexes.clear();
  IdPartStateCombis.clear();
  CombiStateIndexes.clear();
  N=-1;
  J=-1;
  confignum=0;
}

CLinearBasisQN2::CLinearBasisQN2(const CLinearBasisQN2& copybas)
{
  int i;

  Basis=copybas.Basis;
  ParticleConfigs=copybas.ParticleConfigs;
  ConfigStateIndexes=copybas.ConfigStateIndexes;
  IdPartStateCombis=copybas.IdPartStateCombis;
  CombiStateIndexes=copybas.CombiStateIndexes;
  ptypes=copybas.ptypes;
  maxns=copybas.maxns;
  limits=copybas.limits;
  parities=copybas.parities;
  N=copybas.N;
  confignum=copybas.confignum;
  configindex=copybas.configindex;
  if (confignum>0) {
    configurations=new int[confignum*ptypes.size()];
    for (i=0;i<confignum*ptypes.size();i++) {
      configurations[i]=copybas.configurations[i];
    }
  }
  idpartqn=copybas.idpartqn;
  indexlists=copybas.indexlists;
  J=copybas.J;
}

CLinearBasisQN2::~CLinearBasisQN2()
{
  if (confignum>0) { delete [] configurations; };
}

CLinearBasisQN2& CLinearBasisQN2::operator=(const CLinearBasisQN2& assignbas)
{
  int i;

  if (&assignbas!=this) { // no self assignment!
    if (confignum>0) { delete [] configurations; };
    Basis=assignbas.Basis;
    ParticleConfigs=assignbas.ParticleConfigs;
    ConfigStateIndexes=assignbas.ConfigStateIndexes;
    IdPartStateCombis=assignbas.IdPartStateCombis;
    CombiStateIndexes=assignbas.CombiStateIndexes;
    ptypes=assignbas.ptypes;
    maxns=assignbas.maxns;
    limits=assignbas.limits;
    parities=assignbas.parities;
    N=assignbas.N;
    confignum=assignbas.confignum;
    configindex=assignbas.configindex;
    if (confignum>0) {
      configurations=new int[confignum*ptypes.size()];
      for (i=0;i<confignum*ptypes.size();i++) {
	configurations[i]=assignbas.configurations[i];
      }
    }
    idpartqn=assignbas.idpartqn;
    indexlists=assignbas.indexlists;
    J=assignbas.J;    
  }

  return *this;
}

void CLinearBasisQN2::ConstructBasis(int newn,int newj1,int newj2,int newparity,
				    vector<int> newptypes,vector<int> newmaxns,
				    vector<int> newlimits,vector<int> newps)
{
  int i;
  vector<int> newseniorities;

  newseniorities.clear();
  for (i=0;i<newptypes.size();i++) { newseniorities.push_back(-1); };
  ConstructBasis(newn,newj1,newj2,newparity,newptypes,newmaxns,newlimits,
		 newps,newseniorities);
}

void CLinearBasisQN2::ConstructBasis(int newn,int newj1,int newj2,int newparity,
				    vector<int> newptypes,vector<int> newmaxns,
				    vector<int> newlimits,vector<int> newps,
				    vector<int> newseniorities)
{
  int i;
  vector<int> tmp;
  vector<vector<int> > nnsens;

  nnsens.clear();
  for (i=0;i<newseniorities.size();i++) {
    tmp.clear();
    if (newseniorities[i]!=-1) { tmp.push_back(newseniorities[i]); };

    // next line is experimental!
    //    if (newseniorities[i]==-1) { tmp.push_back(newmaxns[i]); };

    nnsens.push_back(tmp);
  }
  ConstructBasis(newn,newj1,newj2,newparity,newptypes,newmaxns,newlimits,
		 newps,nnsens);  
}

void CLinearBasisQN2::ConstructBasis(int newn,int newj1,int newj2,int newparity,
				    vector<int> newptypes,vector<int> newmaxns,
				    vector<int> newlimits,vector<int> newps,
				    vector<vector<int> > newseniorities)
{
  // 'newn'  : particle number of the states
  // 'newj1' and 'newj2' : construct all states with newj1<=jfinal<=newj2
  // 'newparity' : parity of the states, possible values 1,-1 and 0 for both
  // 'newptypes" : the symmetries of the particles
  // 'newmaxns" : maximum particle number of the different particles
  // 'newlimits' : some restictions to the basis
  // 'newps' : the parities of the different particles. 1 or -1
  // 'newseniorities' : allowed seniorities for the id. particle states

  int i,j,k,l,m,n,flag;
  int *configuration,*configarray,*minoccarray,*maxoccarray,*weight;
  int totalconfnum,count,statenum,particlenum;

  //  int nnegmax;

  vector<int> ilist,emptyivec;
  int *js,*limitsptr,*pptr;

  int tmpmaxnu;

  emptyivec.clear();

  N=newn;
  ptypes=newptypes;
  maxns=newmaxns;
  limits=newlimits;
  parities=newps;

  limitsptr=&limits.front();
  pptr=&parities.front();

  Basis.clear();
  ParticleConfigs.clear();
  ConfigStateIndexes.clear();
  IdPartStateCombis.clear();
  CombiStateIndexes.clear();

  statenum=ptypes.size();
  particlenum=N;

  configuration=new int[statenum];
  minoccarray=new int[statenum];
  maxoccarray=new int[statenum];
  weight=new int[statenum];

  for (i=0;i<statenum;i++) {
    weight[i]=1;
    minoccarray[i]=0;
    maxoccarray[i]=maxns[i];
  }

  //  printf("starting to find configurations ..\n");

  count=0;
  totalconfnum=FindConfigurations(particlenum,statenum,minoccarray,
				  maxoccarray,configuration,statenum,
				  0,weight,0,configarray,&count);

  configarray=new int[totalconfnum*statenum];
  count=0;
  totalconfnum=FindConfigurations(particlenum,statenum,minoccarray,
				  maxoccarray,configuration,statenum,
				  0,weight,1,configarray,&count);

  /*
  // print out all the configurations
  for (i=0;i<totalconfnum;i++) {
    printf("[%d] : ",i);
    for (j=0;j<statenum;j++) {
      printf("%d ",configarray[i*statenum+j]);
    }
    printf("\n");
  }
  */

  // Take each configuration and check whether it voilates some
  // restriction of the model. If not, copy it to the global array
  // of all configurations.
  for (j=0;j<2;j++) {
    confignum=0;
    for (i=0;i<totalconfnum;i++) {
      if (LimitsOK(configarray,limitsptr,limits.size()/(ptypes.size()+2),
		   ptypes.size()+2,i)==1) {
	if ((GetParity(configarray+i*statenum,pptr,statenum)==newparity) ||
	    (newparity==0)) {
	  if (j==1) {
	    for (k=0;k<statenum;k++) {
	      configurations[confignum*statenum+k]=configarray[i*statenum+k];
	    }
	  }
	  confignum=confignum+1;
	}
      }
    }
    if (j==0) {
      if (confignum>0) { configurations=new int[confignum*statenum]; };
    }
  }

  // clean up
  delete [] configuration;
  delete [] weight;
  delete [] minoccarray;
  delete [] maxoccarray;
  delete [] configarray;

  /*
  // print out all valid configurations
  for (i=0;i<confignum;i++) {
    printf("[%d] : ",i);
    for (j=0;j<statenum;j++) {
      printf("%d ",configurations[i*statenum+j]);
    }
    printf("\n");
  }
  */

  //  exit(0);

  //  printf("starting to construct id-part. quantum numbers ..\n");

  // construct the identical particle state quantum numbers
  idpartqn.clear();

  for (i=0;i<ptypes.size();i++) {
    idpartqn.push_back(CidParticleStateQN2());

    //    idpartqn[i].ConstructQN(maxns[i],ptypes[i]);
    if (newseniorities[i].size()>0) {
      tmpmaxnu=0;
      for (j=0;j<newseniorities[i].size();j++) {
	if (newseniorities[i][j]>tmpmaxnu) { tmpmaxnu=newseniorities[i][j]; };
      }
      idpartqn[i].ConstructQN(maxns[i],tmpmaxnu,ptypes[i]);
    } else {
      idpartqn[i].ConstructQN(maxns[i],ptypes[i]);
    }
    //    idpartqn[i].Print();
  }


  // Now go through all configurations of particle numbers.
  // For each set of particle numbers, take all possible set of identical
  // particle states. Couple each of theses sets of identical particle
  // states to the final spin and store every possible coupling sheme.

  //  printf("starting to get all possible couplings ..\n");

  indexes=new int[ptypes.size()];
  js=new int[ptypes.size()];
  for (configindex=0;configindex<confignum;configindex++) {
    indexlists.clear();
    flag=1;
    for (j=0;j<ptypes.size();j++) {
      ilist=idpartqn[j].
	GetStatesByQNS(configurations[configindex*statenum+j],
		       newseniorities[j],emptyivec,emptyivec);
      //	GetStatesByQNS(configurations[configindex*statenum+j],
      //		       newseniorities[j],-1,-1);
      //      ilist=idpartqn[j].
      //      	GetStatesByQNS(configurations[configindex*statenum+j],-1,-1,-1);
      indexlists.push_back(ilist);
      if (ilist.size()==0) { flag=0; }; // no combinations this time
    }
    // Now we have ptypes.size() lists of indexes in 'indexlists'.
    // Take every combination  of these states and couple them to each of
    // the final spins in every possible way.

    //      printf("for configuration %d :\n",i);

    if (flag==1) { // if there is at least one combination

      l=Basis.size();
      n=IdPartStateCombis.size();

      for (k=0;k<ptypes.size();k++) { indexes[k]=0; };
      while (indexes[0]<indexlists[0].size()) {
      // Take the current combination and find all possible coupling shemes

      //	printf("combination : ");
      //	for (k=0;k<ptypes.size();k++) { printf("%d ",indexes[k]); };
      //	printf("\n");

	m=Basis.size();
	//	for (J=newj1;J<=newj2;J=J+2) {
	for (J=newj1;J<=newj2;J=J+1) {
	  for (k=0;k<statenum;k++) {
	    // a bug !!!
	    //	    js[k]=idpartqn[k].
	    //	      GetSpin(configurations[configindex*statenum+k],indexes[k]);
	    // a bug !!!
	    js[k]=idpartqn[k].
	      GetSpin(configurations[configindex*statenum+k],
		      indexlists[k][indexes[k]]);
	  }
	  GetAllCouplingShemes(js,ptypes.size(),J);
	}

	if (Basis.size()!=m) { // there are new states for this combination
	  // update arrays
	  for (k=0;k<statenum;k++) {
	    IdPartStateCombis.push_back(indexes[k]);
	  }
	  CombiStateIndexes.push_back(m/(statenum*5));
	  CombiStateIndexes.push_back(Basis.size()/(statenum*5)-1);
	}

	// Go to next combination
	k=ptypes.size()-1;
	indexes[k]=indexes[k]+1;
	while ((indexes[k]>=indexlists[k].size()) && (k>0)) {
	  indexes[k]=0;
	  k=k-1;
	  if (k>=0) { indexes[k]=indexes[k]+1; };
	}
      }

      if (Basis.size()!=l) { // there are new states, update the arrays
	// for the configurations and correspondig indexes
	for (k=0;k<statenum;k++) {
	  ParticleConfigs.push_back(configurations[configindex*statenum+k]);
	}
	ConfigStateIndexes.push_back(l/(statenum*5));
	ConfigStateIndexes.push_back(Basis.size()/(statenum*5)-1);
	ConfigStateIndexes.push_back(n/statenum);
	ConfigStateIndexes.push_back(IdPartStateCombis.size()/statenum-1);
      }

    }

  }

  // clean up
  delete [] indexes;
  delete [] js;
  indexlists.clear();
  idpartqn.clear();
  confignum=0;
  delete [] configurations;
}

void CLinearBasisQN2::GetAllCouplingShemes(int *js,int jnum,int J)
{
  // 'js' is a counter to an array which hold the angular momentums in terms
  // of hbar/2. 'jnum' is the number of these angular momentuns and
  // 'J' is the final spin in terms of hbar/2.

  int *sheme;
  int jrest,jdiff,i,j;

  // sheme[i] x js[i] -> sheme[i+1]

  sheme=new int[jnum+1];

  // construct the initial coupling sheme
  sheme[0]=0;
  for (i=1;i<=jnum;i++) {
    jrest=0;
    for (j=i;j<jnum;j++) {
      jrest=jrest+js[j];
    }
    jdiff=sheme[i-1]-J+js[i-1];
    if (jdiff<=jrest) {
      sheme[i]=sheme[i-1]+js[i-1];
    } else {
      sheme[i]=sheme[i-1]+js[i-1]-(jdiff-jrest);
    }
  }

  // check the initial coupling sheme
  if (!(sheme[jnum]==J)) {
    //      printf("Coupling not possible. Aborting ..\n");
      delete [] sheme;
      return;
  }
  for (i=0;i<jnum;i++) {
    if (!((sheme[i+1]<=sheme[i]+js[i]) &&
	  (abs(sheme[i]-js[i])<=sheme[i+1]) &&
	  ((sheme[i]+js[i])%2==sheme[i+1]%2))) {
      //      printf("Coupling not possible. Aborting ..\n");
      delete [] sheme;
      return;
    }
  }

  //  PrintSheme(sheme,js,jnum);
  StoreStateQN(sheme,js,jnum);

  // now go through all possible coupling shemes
  while (GetNextSheme(js,sheme,jnum,J)==1) {
    //    PrintSheme(sheme,js,jnum);
    StoreStateQN(sheme,js,jnum);
  }

  //  delete [] js;
  delete [] sheme;
}

int CLinearBasisQN2::GetNextSheme(int* js,int* sheme,int jnum,int J)
{
  // Take the actual coupling sheme given by the arguments and go to
  // the next one. The array 'sheme' will be changed. In case of success
  // this functions returns '1' and '0' otherwise. In case of '0' there
  // is no next coupling sheme.

  int i,j,k,res,jrest,jdiff;

  i=jnum-1;
  res=0;
  while (i>0) {
    sheme[i]=sheme[i]-2;
    // check if coupling to the right side is still ok
    if (!((sheme[i+1]<=sheme[i]+js[i]) &&
	  (abs(sheme[i]-js[i])<=sheme[i+1]))) {
      // coupling not ok, undo and go to lower index
      sheme[i]=sheme[i]+2;
      i=i-1;
    } else { // coupling to right side ok, check left side
      if (!((sheme[i]<=sheme[i-1]+js[i-1]) &&
	    (abs(sheme[i-1]-js[i-1])<=sheme[i]))) {
	i=i-1;
      } else { // coupling ok
	for (j=i+1;j<=jnum;j++) {
	  jrest=0;
	  for (k=j;k<jnum;k++) {
	    jrest=jrest+js[k];
	  }
	  jdiff=sheme[j-1]-J+js[j-1];
	  if (jdiff<=jrest) {
	    sheme[j]=sheme[j-1]+js[j-1];
	  } else {
	    sheme[j]=sheme[j-1]+js[j-1]-(jdiff-jrest);
	  }
	  // check if coupling is OK
	  if ((sheme[j]<=sheme[j-1]+js[j-1]) &&
	      (abs(sheme[j-1]-js[j-1])<=sheme[j]) &&
	      ((sheme[j-1]+js[j-1])%2==sheme[j]%2) &&
	      (sheme[jnum]==J)) {
	    res=1; // take this sheme and
	    i=0;   // exit while loop
	  }
	}
      }
    }
  }

  return res;
}

void CLinearBasisQN2::StoreStateQN(int* sheme,int* js,int jnum)
{
  // Add a new state to the vector 'Basis'. The neccessary informations
  // are taken from the arguments as well as from
  // the gobal data in the private section.
  int i;
  int N,nu,delta;

  for (i=0;i<jnum;i++) {
    N=configurations[configindex*jnum+i];
    nu=idpartqn[i].GetNu(N,indexlists[i][indexes[i]]);
    delta=idpartqn[i].GetDelta(N,indexlists[i][indexes[i]]);
    Basis.push_back(N);
    Basis.push_back(nu);
    Basis.push_back(delta);
    Basis.push_back(js[i]);
    Basis.push_back(sheme[i+1]);
  }
}

/*
void CLinearBasisQN::FilterSeniority(vector<int> seniorities)
{
  // Filter out states with some fixed seniorities. All other states
  // will be erased from the basis.
  // Arguments: 'seniorities' is an array of ptypes.size() integers.
  // For each particle type it gives the seniority of the states which
  // should NOT be deleted. A '-1' stands for 'all seniorities'.
  vector<int> newbas;
  int i,j,statenum,statesize,flag;

  if (seniorities.size()!=ptypes.size()) {
    printf("void CLinearBasisQN::FilterSeniority(vector<int>) :\n");
    printf("Number of given seniorities do not match the number of\n");
    printf("particle types. Aborting..\n");
    exit(0);
  }

  newbas.clear();
  statenum=GetStateNum();
  statesize=ptypes.size()*5;
  for (i=0;i<statenum;i++) {
    flag=1;
    for (j=0;j<seniorities.size();j++) {
      if ((Basis[i*statesize+j*5+1]!=seniorities[j]) && (seniorities[j]!=-1)) {
	flag=0;
      }
    }
    if (flag==1) { // take that state
      for (j=0;j<statesize;j++) {
	newbas.push_back(Basis[i*statesize+j]);
      }
    }
  }
  Basis=newbas;
}
*/

void CLinearBasisQN2::Print()
{
  // print out the quantum numbers of all states
  int i,j;

  if (ptypes.size()==0) { return; };
  for (i=0;i<Basis.size()/(ptypes.size()*5);i++) {
    for (j=0;j<ptypes.size();j++) {
      if (j==0) {
	printf("|%d> = ",i);
      } else {
	printf(" x ");
      }
      printf("|%d,%d,%d,%d/2>(%d/2)",Basis[i*ptypes.size()*5+j*5+0],
	     Basis[i*ptypes.size()*5+j*5+1],Basis[i*ptypes.size()*5+j*5+2],
	     Basis[i*ptypes.size()*5+j*5+3],Basis[i*ptypes.size()*5+j*5+4]);
    }
    printf("\n");
  }
}

int CLinearBasisQN2::GetStateNum()
{
  int res;
  if (ptypes.size()<=0) {
    printf("int CLinearBasisQN::GetStateNum()\n");
    printf("Number of different particles <= 0. Aborting..\n");
    exit(0);
  }
  res=Basis.size()/(ptypes.size()*5);
  return res;
}

int CLinearBasisQN2::GetConfNum()
{
  // Return the number of configurations
  int res;
  if (ptypes.size()<=0) {
    printf("int CLinearBasisQN::GetStateNum()\n");
    printf("Number of different particles <= 0. Aborting..\n");
    exit(0);
  }
  res=ParticleConfigs.size()/ptypes.size();
  return res;
}

void CLinearBasisQN2::Clear()
{
  // clean up everything
  Basis.clear();
  ptypes.clear();
  maxns.clear();
  limits.clear();
  ParticleConfigs.clear();
  ConfigStateIndexes.clear();
  IdPartStateCombis.clear();
  CombiStateIndexes.clear();
  N=-1;
  J=-1;
  if (confignum>0) { delete [] configurations; };
  confignum=0;
}

void CLinearBasisQN2::PrintConfigs()
{
  int i,j;

  if (ptypes.size()==0) { return; };
  for (i=0;i<ParticleConfigs.size()/ptypes.size();i++) {
    printf("[%d] : ",i);
    for (j=0;j<ptypes.size();j++) {
      printf("%d ",ParticleConfigs[i*ptypes.size()+j]);
    }
    printf(" for states %d,..,%d and combis %d,..,%d\n",
	   ConfigStateIndexes[i*4],ConfigStateIndexes[i*4+1],
	   ConfigStateIndexes[i*4+2],ConfigStateIndexes[i*4+3]);
  }
}

void CLinearBasisQN2::PrintCombis()
{
  int i,j;

  if (ptypes.size()==0) { return; };
  for (i=0;i<IdPartStateCombis.size()/ptypes.size();i++) {
    printf("(%d) : ",i);
    for (j=0;j<ptypes.size();j++) {
      printf("%d ",IdPartStateCombis[i*ptypes.size()+j]);
    }
    printf(" belongs to states : %d,..,%d\n",
	   CombiStateIndexes[i*2],CombiStateIndexes[i*2+1]);
  }
}

int* CLinearBasisQN2::GetPtr(int stateindex)
{
  // return a pointer to the state with index 'stateindex'

  int *res;

  res=&Basis.front()+stateindex*ptypes.size()*5;

  return res;
}

void CLinearBasisQN2::GetJArray(vector<int>* jarray)
{
  int* ptr;
  int statenum,ptnum,i;

  ptr=&Basis.front();
  statenum=GetStateNum();
  ptnum=ptypes.size();

  jarray->clear();

  for (i=0;i<statenum;i++) {
    jarray->push_back(ptr[(i+1)*ptnum*5-1]);
  }
}

void CLinearBasisQN2::SelectStates(int* stateptr,int optype,vector<int>* statelist)
{
  int *neededconf,*basptr;
  int i,j,ptnum,ci,flag,seniority1,seniority2;
  int j1,j1p,j2,j2p,J,Jp,opj;  // for some angular momenta

  if (optype==0) {
    printf("optype=0 !?!? Aborting ..\n");
    exit(0);
  }

  statelist->clear();
  ptnum=ptypes.size();

  //  printf("state : ");
  //  for (i=0;i<ptypes.size()*5;i++) { printf("%d ",stateptr[i]); };
  //  printf("\n");

  neededconf=new int[ptypes.size()];

  basptr=&Basis.front();

  for (i=0;i<ptypes.size();i++) {
    neededconf[i]=stateptr[5*i];
  }

  if (optype>0) {
    neededconf[abs(optype)-1]=neededconf[abs(optype)-1]-1;
  } else {
    neededconf[abs(optype)-1]=neededconf[abs(optype)-1]+1;
  }

  //  printf("needed conf : ");
  //  for (i=0;i<ptypes.size();i++) { printf("%d ",neededconf[i]); };
  //  printf("\n");

  // Now get the index of the configuration which is stored in 'neededconf'.
  // -1 means that this configuration does not exist.
  ci=GetConfIndex(neededconf);

  //  printf("index of needed conf : %d\n",ci);

  if (ci>=0) {
    //    printf("interesting states due to particle numbers : %d,..,%d\n",
    //	   ConfigStateIndexes[ci*4],ConfigStateIndexes[ci*4+1]);

    // Check the selection rules with all these states and store the indexes
    // of the remaining states in an array.
    for (i=ConfigStateIndexes[ci*4];i<=ConfigStateIndexes[ci*4+1];i++) {
      //      printf("checking state %d : ",i);
      //      for (j=0;j<ptnum*5;j++) { printf("%d ",basptr[i*ptnum*5+j]); };
      //      printf("\n");
      j=0;
      flag=0;
      while ((j<(abs(optype)-1)*5) && (flag==0)) {
	if (stateptr[j]!=basptr[i*ptnum*5+j]) {
	  flag=1;
	}
	j=j+1;
      }

      if (flag==0) {
	//	printf("first delta-selection-rule is ok for state : %d\n",i);
       	j=abs(optype)*5;
	flag=0;
	while ((j<ptnum*5) && (flag==0)) {
	  // Don't check the intermediate coupling. That we do in the netxt step
	  if ((j+1)%5!=0) {
	    if (stateptr[j]!=basptr[i*ptnum*5+j]) {
	      flag=1;
	    }
	  }
	  j=j+1;
	}
	if (flag==0) {
	  //	  printf("last delta-selection-rule is ok for state : %d\n",i);

	  // During the process of the successiv application of
	  // <a1,j1,a2,j2,J||T(1)^k||a1',j1',a2',j2',J'> = (-)^(j1+j2+J'+k) *
	  // sqrt((2J+1)(2J'+1)) * <a1,j1||T(1)^k||a1',j1'> *
	  // wigner6j(j1,J,j2,J',j1',k) * delta(a2,a2') * delta(j2,j2')
	  // we have to check for each application of this rule four
	  // angular momentum selection rules according to the wigner6j symbol.
	  j=ptnum*5-1; // index to the last (to the final) angular momentum
	  flag=0;
	  while ((j>abs(optype)*5) && (flag==0)) {
	    j1=stateptr[j-5];
	    j2=stateptr[j-1];
	    J=stateptr[j];
	    j1p=basptr[i*ptnum*5+j-5];
	    j2p=basptr[i*ptnum*5+j-1];
	    Jp=basptr[i*ptnum*5+j];
	    opj=ptypes[abs(optype)-1]-1;
	    //	    printf("%d %d %d %d %d %d\n",j1,j2,J,opj,j1p,j2p,Jp);
	    // Now check the four rules for wigner6j(j1,J,j2,Jp,j1p,opj)
	    // These rules are: (j1,J,j2),(Jp,j1p,j2),(j1,j1p,opj),(Jp,J,opj)
	    if (!((abs(j1-J)<=j2) && (j2<=j1+J))) { flag=1; };
	    if (!((abs(Jp-j1p)<=j2) && (j2<=Jp+j1p))) { flag=1; };
	    if (!((abs(j1-j1p)<=opj) && (opj<=j1+j1p))) { flag=1; };
	    if (!((abs(Jp-J)<=opj) && (opj<=Jp+J))) { flag=1; };
	    j=j-5; // check next wigner6j symbol
	  }
	  if (flag==0) {
	    //	    printf("wigner6j-selection rules are ok for state : %d\n",i);
	    // Check seniority selection rule
	    seniority1=basptr[i*ptnum*5+(abs(optype)-1)*5+1];
	    seniority2=stateptr[(abs(optype)-1)*5+1];
	    if (abs(seniority1-seniority2)==1) {
	      // printf("seniority-selection-rule is ok for state : %d\n",i);
	      // Check angular momentum selection rule

	      /*
	      j1=basptr[i*ptnum*5+(abs(optype)-1)*5+3];
	      j2=stateptr[(abs(optype)-1)*5+3];
	      opj=ptypes[abs(optype)-1]-1;
	      if ((j1<=j2+opj) && (j1>=abs(j2-opj))) {
		// printf("angular-momentum selection rule is ok for state %d\n",i);
		statelist->push_back(i);
	      }
	      */

	      flag=0;
	      if (abs(optype)>1) {
		j1=basptr[i*ptnum*5+(abs(optype)-1)*5-1];
		j2=basptr[i*ptnum*5+(abs(optype)-1)*5+3];
		J=basptr[i*ptnum*5+(abs(optype)-1)*5+4];
		j1p=stateptr[(abs(optype)-1)*5-1];
		j2p=stateptr[(abs(optype)-1)*5+3];
		Jp=stateptr[(abs(optype)-1)*5+4];
		opj=ptypes[abs(optype)-1]-1;
		// Check wigner6j(j2,J,j1,Jp,j2p,opj) selection rules.
		// These are: (j2,J,j1),(Jp,j2p,j1),(j2,j2p,opj),(Jp,J,opj)
		if (!((abs(j2-J)<=j1) && (j1<=j2+J))) { flag=1; };
		if (!((abs(Jp-j2p)<=j1) && (j1<=Jp+j2p))) { flag=1; };
		if (!((abs(j2-j2p)<=opj) && (opj<=j2+j2p))) { flag=1; };
		if (!((abs(Jp-J)<=opj) && (opj<=Jp+J))) { flag=1; };
	      }
	      if (abs(optype)==1) {
		j1=basptr[i*ptnum*5+(abs(optype)-1)*5+3];
		j2=stateptr[(abs(optype)-1)*5+3];
		opj=ptypes[abs(optype)-1]-1;
		if (!((j1<=j2+opj) && (j1>=abs(j2-opj)))) { flag=1; };
	      }
	      if (flag==0) {
		// printf("angular-momentum selection rule is ok for state %d\n",i);
		statelist->push_back(i);
	      }
	    }
	  }
	}
      }
    }
  }

  //  if (statelist->size()>0) {
  //    printf("statelist : ");
  //    for (i=0;i<statelist->size();i++) { printf("%d ",(*statelist)[i]); };
  //    printf("\n");
  //  }

  /*
  for (j=0;j<statelist->size();j++) {
    printf("M: <");
    for (i=0;i<ptnum*5;i++) { printf("%d ",stateptr[i]); };
    printf("|| %d || ",optype);
    for (i=0;i<ptnum*5;i++) { printf("%d ",basptr[(*statelist)[j]*ptnum*5+i]); };
    printf(">\n");
  }
  */

  delete [] neededconf;
}

int CLinearBasisQN2::GetConfIndex(int* searchconf)
{
  int i,j,k,res,flag,confnum,ptnum,*confptr;

  ptnum=ptypes.size();

  // If there a negative particle number then return a 'not found'
  flag=0;
  for (i=0;i<ptnum;i++) {
    if (searchconf[i]<0) { flag=1; };
  }
  if (flag==1) { return -1; }; // -1 means 'not found'

  flag=0;
  confptr=&ParticleConfigs.front();
  confnum=GetConfNum();
  res=-1;
  i=0;
  while (flag==0) {
    if (i<confnum) {
      j=0;
      k=1;
      for (j=0;j<ptnum;j++) {
	if (searchconf[j]!=confptr[i*ptnum+j]) {
	  k=0;
	  j=ptnum; // exit for-loop
	}
      }
      if (k==1) {
	flag=1;
	res=i;
      }
    } else {
      flag=1;
    }
    i=i+1;
  }

  return res;
}

////////////////////////////////////////////////
/////    ** private functions **
////////////////////////////////////////////////

void CLinearBasisQN2::PrintSheme(int* sheme,int* js,int jnum)
{
  // this function has only debugging purpose
  int i;
  printf("%d/2 x %d/2 ",sheme[0],js[0]);
  for (i=1;i<jnum;i++) {
    printf("-> %d/2 x %d/2 ",sheme[i],js[i]);
  }
  printf("-> %d/2\n",sheme[jnum]);
}

int CLinearBasisQN2::FindConfigurations(int particlenum,int statenum,
				       int* minoccupation,int* maxoccupation,
				       int* configuration,int initialstatenum,
				       int writeconfig,int* weight,int storeflag,
				       int* configarray,int* count)
{
  // This is a recursiv algorithm
  // It goes through all possible occupations of the first state
  // and then calls itself with the remaining particles and states

  // parameters :
  // particlenum     = number of identical particles
  // statenum        = number of states the particles can occupy
  // minoccupation[] = minimum number of particles in each state
  // maxoccupation[] = maximum number of partciles in each state
  // configuration   = pointer to an integer array to store the configuration
  //                   The size must be grater ot equal initialstatenum
  // initialstatenum = In the initial call is statenum==initialstatenum.
  //                   Initialstatenum is needed internally
  // writeconfig     = 1: all configurations are printed out
  //                   0: the configurations will not be printed out
  // weight          = Array with a weight for each state. Thats difficult
  // to explain : For the finding of simple distribution of particles over
  // a number of states you have to fill this array with '1'.
  // storeflag       = 1: all configurations will be stored in configarray[]
  //                      this array has to be already allocated
  //                   0: the configuration will not be stored
  // configarray     = pointer to an array of appropriate size to
  //                   receive all configurations (if storeflag==1)
  // count           = pointer to an integer which must be 0 at initial call
  // The return value is the number of possible configurations

  int confignum=0; // count the configurations found
  int minfirststateoccupation; // the minimum and maximum occupation
  int maxfirststateoccupation; // number possible for the first state
  int occupation; // occupation number of the first state
  int i;
  int summinoccupation,summaxoccupation;

  // break conditions are :
  // statenum = 0
  // \sum minoccupation[i] > partcilenum
  // \sum maxoccupation[i] < particlenum
  summinoccupation=0;
  summaxoccupation=0;
  for (i=0;i<statenum;i++) {
    summinoccupation=summinoccupation+minoccupation[i];
    summaxoccupation=summaxoccupation+maxoccupation[i];
  }
  if ((statenum>0) && (summinoccupation<=particlenum) &&
      (summaxoccupation>=particlenum)) {

    // calculate the possible occupation numbers for the first state
    minfirststateoccupation=minoccupation[statenum-1];
    // maxfirststateoccupation=min(particlenum,maxoccupation)
    maxfirststateoccupation=particlenum;
    if (maxoccupation[statenum-1]<particlenum) {
      maxfirststateoccupation=maxoccupation[statenum-1];
    }

    if (statenum==1) {
      confignum=confignum+1; // only one possible occupation
      configuration[statenum-1]=particlenum; // store the occupation

      if (writeconfig==1) {
	// write out the configuration
	cout << "configuration " << (*count) << " : ";
	for (i=0;i<initialstatenum;i++) {
	  cout << " " << configuration[i];
	}
	cout << endl;
      }
      if (storeflag==1) {
	for (i=0;i<initialstatenum;i++) {
	  configarray[(*count)*initialstatenum+i]=configuration[i];
	}
      }
      (*count)=(*count)+1;
    } else {
      // go through all possible occupations of the first state
      for (occupation=minfirststateoccupation;
	   occupation<=maxfirststateoccupation;occupation++) {
	configuration[statenum-1]=occupation; // store the occupation
	// recursiv call for the remaining partciles and states
	confignum=confignum+
	  FindConfigurations(particlenum-occupation*weight[statenum-1],
			     statenum-1,minoccupation,
			     maxoccupation,configuration,
			     initialstatenum,writeconfig,weight,
			     storeflag,configarray,count);
      }
    }
  }

  //  cout << particlenum << " " << statenum << " ";
  //  cout << minoccupation << " " << maxoccupation << " ";
  //  cout << confignum << endl;
  return confignum; // give back the number of configurations found
}

int CLinearBasisQN2::LimitsOK(int* confs,int* lims,int limnum,
			     int limlength,int iconf)
{
  // 'lims'      : pointer to the array with the restrictions
  // 'limnum'    : number of the restrictions
  // 'limlength' : number of integers of one restrictions
  // 'iconf'     : index of the configuration to check

  int i,j,k,res,ptnum;

  ptnum=ptypes.size();

  res=1;
  for (i=0;i<limnum;i++) {
    k=0;
    for (j=0;j<limlength-2;j++) {
      k=k+lims[i*limlength+j]*confs[iconf*ptnum+j];
    }
    if (lims[(i+1)*limlength-2]==1) { // check for equal
      if (!(lims[(i+1)*limlength-1]==k)) { res=0; };
    }
    if (lims[(i+1)*limlength-2]==2) { // check for less or equal
      if (!(k<=lims[(i+1)*limlength-1])) { res=0; };
    }
  }

  return res;
}

int CLinearBasisQN2::GetParity(int* conf,int* ps,int ptnum)
{
  int res,i;
  res=1;
  for (i=0;i<ptnum;i++) {
    if (conf[i]%2==1) {
      if (ps[i]==-1) {
	res=-res;
      }
    }
  }
  return res;
}
