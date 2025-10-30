////////////////////////////////////////////////////////////////////

#include "./cdrecouple.h"

CDRecouple::CDRecouple()
{
  int i;

  ClearAll();

  FCT=new double[FCT_Size];
  FCT[0]=0;
  FCT[1]=1;
  for (i=2;i<FCT_Size;i++) {
    FCT[i]=FCT[i-1]*(double)(i-1)/10.0;
  }
}

CDRecouple::CDRecouple(const CDRecouple& copyterms)
{
  int i;

  FCT=new double[FCT_Size];
  for (i=0;i<FCT_Size;i++) {
    FCT[i]=copyterms.FCT[i];
  }

  ParticleTypes=copyterms.ParticleTypes;
  Termsddtt=copyterms.Termsddtt;
  Termsdt=copyterms.Termsdt;
}

CDRecouple::~CDRecouple()
{
  delete [] FCT;
}

CDRecouple& CDRecouple::operator=(const CDRecouple& assigncdr)
{
  int i;
  if (&assigncdr!=this) {
    delete [] FCT;
    FCT=new double[FCT_Size];
    for (i=0;i<FCT_Size;i++) {
      FCT[i]=assigncdr.FCT[i];
    }
    ParticleTypes=assigncdr.ParticleTypes;
    Termsddtt=assigncdr.Termsddtt;
    Termsdt=assigncdr.Termsdt;
  }
  return *this;
}

void CDRecouple::Add_onebody(int n1,int n2,double factor)
{
  // Add a new term of the form (c+_{n1} x c~_{n2})^0

  vector<int> term;
  map<vector<int>,double>::iterator iter;

  term.clear();
  term.push_back(n1);
  term.push_back(n2);

  if ((iter=Termsdt.find(term))==Termsdt.end()) {
    Termsdt[term]=factor;
  } else {
    Termsdt[term]=Termsdt[term]+factor;
  }
}

void CDRecouple::Add_ddtt_to_ddtt(int n1,int n2,int n3,int n4,
				  int k,double factor)
{
  // Add a new term of the form
  // 'factor'*[ (c+_{n1} x c+_{n2})^k x (c~_{n3} x c~_{n4})^k ]^0.

  vector<int> term;
  map<vector<int>,double>::iterator iter;

  term.clear();
  term.push_back(n1);
  term.push_back(n2);
  term.push_back(n3);
  term.push_back(n4);
  term.push_back(k);
  if ((iter=Termsddtt.find(term))==Termsddtt.end()) {
    Termsddtt[term]=factor;
  } else {
    Termsddtt[term]=Termsddtt[term]+factor;
  }
}

void CDRecouple::Add_dtdt_to_ddtt(int n1,int n2,int n3,int n4,
				  int k,double factor)
{
  // Add a new term of the form
  // 'factor'*[ (c+_{n1} x c~_{n2})^k x (c+_{n3} x c~_{n4})^k ]^0.
  // This term will be recoupled to the form
  // 'factor'*(
  // \sum_J F(J) * [ (c+_{n1} x c+_{n3})^J x (c~_{n2} x c~_{n4})^J ]^0
  // + F1 * one-body )
  // and added to the already stored ones.
  // The values of n1 to n4 correspond to the particle types stored
  // in 'ParticleTypes'. The first particle has n1=1.
  // J has to be given in units of hbar/2.

  double coef;
  int j1,j2,j3,j4;
  int J,phase,jstep;
  vector<int> term;
  map<vector<int>,double>::iterator iter;

  j1=ParticleTypes[n1-1]-1;
  j2=ParticleTypes[n2-1]-1;
  j3=ParticleTypes[n3-1]-1;
  j4=ParticleTypes[n4-1]-1;

  //  printf("[(c+_{%d,%d/2}xc~_{%d,%d/2})^(%d) x ",n1,j1,n2,j2,k);
  //  printf("(c+_{%d,%d/2}xc~_{%d,%d/2})^(%d)]^(0) = \n",n3,j3,n4,j4,k);

  // We start with the onebody term:
  if ((n1==n4) && (n2==n3)) {
    coef=sqrt((double)(k+1)/(j1+1));
    if ((k+j1+j3)%4==0) { phase=1; } else { phase=-1; };
    if (j2%2==1) { phase=-phase; };
    coef=coef*phase;
    //    printf("%f * %f * \n",factor,coef);
    //    printf("[c+_{%d,%d/2} x c~_{%d,%d/2}]^(0)\n",n1,j1,n4,j4);
    // add term to the map 'Termsdt'
    term.clear();
    term.push_back(n1);
    term.push_back(n4);
    if ((iter=Termsdt.find(term))==Termsdt.end()) {
      Termsdt[term]=factor*coef;
    } else {
      Termsdt[term]=Termsdt[term]+factor*coef;
    }
  }

  // No do the two body terms
  if (n1==n3) { jstep=4; } else { jstep=2; };
  for (J=abs(j1-j3);J<=j1+j3;J=J+jstep) {
    if ((J<=j2+j4) && (J>=abs(j2-j4))) {
      coef=sqrt((double)(k+1)*(J+1))*wigner6j(j4,j2,J,j1,j3,k);
      if ((j2+j3+k+J)%4==0) { phase=1; } else { phase=-1; };
      if (n2==n3) { if (j2%2==1) { phase=-phase; }; };
      coef=coef*phase;
      //    printf("%f * %f * \n",factor,coef);
      //    printf("[(c+_{%d,%d/2}xc+_{%d,%d/2})^(%d) x ",n1,j1,n3,j3,J);
      //    printf("(c~_{%d,%d/2}xc~_{%d,%d/2})^(%d)]^(0)\n",n2,j2,n4,j4,J);
      term.clear();
      term.push_back(n1);
      term.push_back(n3);
      term.push_back(n2);
      term.push_back(n4);
      term.push_back(J);
      if ((iter=Termsddtt.find(term))==Termsddtt.end()) {
	Termsddtt[term]=factor*coef;
      } else {
	Termsddtt[term]=Termsddtt[term]+factor*coef;
      }
    }
  }
}

void CDRecouple::SetParticleTypes(vector<int> ptypes)
{
  ParticleTypes=ptypes;
}

void CDRecouple::Print()
{
  vector<int> key;
  double value;
  map<vector<int>,double>::const_iterator iter;
  int n1,n2,n3,n4,j1,j2,j3,j4,J;
  int termcount;

  termcount=0;
  for (iter=Termsdt.begin();iter!=Termsdt.end();iter++) {
    key=iter->first;
    value=iter->second;
    n1=key[0];
    n2=key[1];
    j1=ParticleTypes[n1-1]-1;
    j2=ParticleTypes[n2-1]-1;
    if (termcount>0) { printf(" +\n"); };
    printf("%f * ",value);
    printf("[c+_{%d,%d/2} x c~_{%d,%d/2}]^(0)",n1,j1,n2,j2);
    termcount=termcount+1;
  }
  for (iter=Termsddtt.begin();iter!=Termsddtt.end();iter++) {
    key=iter->first;
    value=iter->second;
    n1=key[0];
    n2=key[1];
    n3=key[2];
    n4=key[3];
    J=key[4];
    j1=ParticleTypes[n1-1]-1;
    j2=ParticleTypes[n2-1]-1;
    j3=ParticleTypes[n3-1]-1;
    j4=ParticleTypes[n4-1]-1;
    if (termcount>0) { printf(" +\n"); };
    printf("%f * ",value);
    printf("[(c+_{%d,%d/2}xc+_{%d,%d/2})^%d x ",n1,j1,n2,j2,J);
    printf("(c~_{%d,%d/2}xc~_{%d,%d/2})^%d]^0",n3,j3,n4,j4,J);
    termcount=termcount+1;
  }
  printf("\n");
}

void CDRecouple::ClearAll()
{
  ParticleTypes.clear();
  Termsddtt.clear();
  Termsdt.clear();
}

int CDRecouple::GetOneBodyNum()
{
  return Termsdt.size();
}

int CDRecouple::GetTwoBodyNum()
{
  return Termsddtt.size();
}

void CDRecouple::GetOneBodyTerms(vector<int>* terms,vector<double>* coefs)
{
  map<vector<int>,double>::iterator iter;
  vector<int> term;

  terms->clear();
  coefs->clear();

  for (iter=Termsdt.begin();iter!=Termsdt.end();iter++) {
    term=iter->first;
    terms->push_back(term[0]);
    terms->push_back(term[1]);
    coefs->push_back(iter->second);
  }
}

void CDRecouple::GetTwoBodyTerms(vector<int>* terms,vector<double>* coefs)
{
  map<vector<int>,double>::iterator iter;
  vector<int> term;

  terms->clear();
  coefs->clear();

  for (iter=Termsddtt.begin();iter!=Termsddtt.end();iter++) {
    term=iter->first;
    terms->push_back(term[0]);
    terms->push_back(term[1]);
    terms->push_back(term[2]);
    terms->push_back(term[3]);
    terms->push_back(term[4]);
    coefs->push_back(iter->second);
  }
}

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

double CDRecouple::wigner6jhelp(int J1,int J2,int J3)
{
  // all in units of hbar/2
  // this is only a sub-function which is needed by the wigner6j-function
  int IW1,IW2,IW3,IW4;
  double FDELTA;
  IW1=(int)(J1+J2-J3)/2+1;
  IW2=(int)(J1-J2+J3)/2+1;
  IW3=(int)(-J1+J2+J3)/2+1;
  IW4=(int)(J1+J2+J3+2)/2+1;
  if ((IW1>=FCT_Size) || (IW2>=FCT_Size) || (IW3>=FCT_Size) ||
      (IW4>=FCT_Size)) {
    cout << "double CDidParticleModel::wigner6jhelp(int,int,int) :" << endl;
    cout << "Index out of array range, you have to ";
    cout << "enlarge FCT[] !" << endl;
    exit(0);
  }
  FDELTA=sqrt(FCT[IW1]*FCT[IW2]*FCT[IW3]/FCT[IW4]);
  return FDELTA/sqrt((double)10.0);
}

double CDRecouple::wigner6j(int J1,int J2,int J3,int L1,int L2,int L3)
{
  // / J1 J2 J3 \
  // \ L1 L2 L3 /
  // all in units of hbar/2
  //  int i;

  double OMEGA,PH;
  int IWMIN,IWMAX,IW,IW1,IW2,IW3,IW4,IW5,IW6,IW7,IW8;
  double CC=0;
  if (J1+J2-J3>=0) {
    if (abs(J1-J2)-J3<=0) {
      if ((J1+J2+J3-2*((int)((J1+J2+J3)/2)))==0) {
        if (J1+L2-L3>=0) {
          if (abs(J1-L2)-L3<=0) {
            if ((J1+L2+L3-2*((int)(J1+L2+L3)/2))==0) {
              if (L1+J2-L3>=0) {
                if (abs(L1-J2)-L3<=0) {
                  if ((L1+J2+L3-2*((int)(L1+J2+L3)/2))==0) {
                    if (L1+L2-J3>=0) {
                      if (abs(L1-L2)-J3<=0) {
                        if ((L1+L2+J3-2*((int)(L1+L2+J3)/2))==0) {
                          OMEGA=0;
                          if (J3==0) {
                            return pow(-1,(double)(J1+L2+L3)/2)/
                              sqrt(((double)(J1)+1)*((double)(L2)+1));
                          }
                          if (L3==0) {
                            return pow(-1,(double)(J1+J2+J3)/2)/
                              sqrt(((double)(J1)+1)*((double)(J2)+1));
                          }
                          IWMIN=J1+J2+J3;
                          if (IWMIN-J1-L2-L3<0) {
                            IWMIN=J1+L2+L3;
                          }
                          if (IWMIN-L1-J2-L3<0) {
                            IWMIN=L1+J2+L3;
                          }
                          if (IWMIN-L1-L2-J3<0) {
                            IWMIN=L1+L2+J3;
                          }
                          IWMAX=J1+J2+L1+L2;
                          if (IWMAX-J2-J3-L2-L3>0) {
                            IWMAX=J2+J3+L2+L3;
                          }
                          if (IWMAX-J1-J3-L1-L3>0) {
                            IWMAX=J1+J3+L1+L3;
                          }
                          if (IWMIN-IWMAX<=0) {
                            for (IW=IWMIN;IW<=IWMAX;IW=IW+2) {
                              IW1=(int)IW/2+2;
                              IW2=(int)(IW-J1-J2-J3)/2+1;
                              IW3=(int)(IW-J1-L2-L3)/2+1;
                              IW4=(int)(IW-L1-J2-L3)/2+1;
                              IW5=(IW-L1-L2-J3)/2+1;
                              IW6=(J1+J2+L1+L2-IW)/2+1;
                              IW7=(int)(J1+J3+L1+L3-IW)/2+1;
                              IW8=(int)(J2+J3+L2+L3-IW)/2+1;
                              if (IW-4*(int)(IW/4)==0) {
                                PH=1;
                              } else {
                                PH=-1;
                              }
			      if ((IW1>=FCT_Size) || (IW2>=FCT_Size) ||
				  (IW3>=FCT_Size) || (IW4>=FCT_Size) ||
				  (IW5>=FCT_Size) || (IW6>=FCT_Size) ||
				  (IW7>=FCT_Size) || (IW8>=FCT_Size)) {
				cout << "double wigner6";
				cout << "(int,int,int,int,int,int) :" << endl;
				cout << "Index out of array range, ";
				cout << "you have to enlarge FCT[] !" << endl;
				exit(0);
			      }
                              OMEGA=OMEGA+PH*FCT[IW1]/FCT[IW2]/FCT[IW3]/
                                FCT[IW4]/FCT[IW5]/FCT[IW6]/FCT[IW7]/FCT[IW8];
                            }
                            CC=OMEGA*wigner6jhelp(J1,J2,J3)*
                              wigner6jhelp(J1,L2,L3)*
                              wigner6jhelp(L1,J2,L3)*
                              wigner6jhelp(L1,L2,J3);
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  return CC*10.0;  
}

////////////////////////////////////////////////////////////////////
