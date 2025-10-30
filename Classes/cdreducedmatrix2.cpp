//////////////////////////////////////////////////////////////////////

#include "./cdreducedmatrix2.h"

CDReducedMatrix2::CDReducedMatrix2()
{
  int i;

  NonZeros=0;
  NCols=0;
  NRows=0;

  IMatrix=new int[1];
  IMatrix[0]=0;

  FCT=new double[FCT_Size];
  FCT[0]=0;
  FCT[1]=1;
  for (i=2;i<FCT_Size;i++) {
    FCT[i]=FCT[i-1]*(double)(i-1)/10.0;
  }
}

CDReducedMatrix2::CDReducedMatrix2(const CDReducedMatrix2& copymat)
{
  int i;

  NonZeros=copymat.NonZeros;
  NRows=copymat.NRows;
  NCols=copymat.NCols;
  IMatrix=new int[NRows+1];
  for (i=0;i<NRows+1;i++) {
    IMatrix[i]=copymat.IMatrix[i];
  }
  if (NonZeros>0) {
    Matrix=new double[NonZeros];
    JMatrix=new int[NonZeros];
    for (i=0;i<NonZeros;i++) {
      Matrix[i]=copymat.Matrix[i];
      JMatrix[i]=copymat.JMatrix[i];
    }
  }

  FCT=new double[FCT_Size];
  for (i=0;i<FCT_Size;i++) { FCT[i]=copymat.FCT[i]; };
}

CDReducedMatrix2::~CDReducedMatrix2()
{
  ClearMatrix();
  delete [] IMatrix;
  if (FCT_Size>0) { delete [] FCT; };
}

CDReducedMatrix2& CDReducedMatrix2::operator=(const CDReducedMatrix2& eqmat)
{
  int i;

  if (&eqmat!=this) { // no self assignment!
    (*this).ClearMatrix();
    (*this).NonZeros=eqmat.NonZeros;
    (*this).NRows=eqmat.NRows;
    (*this).NCols=eqmat.NCols;
    (*this).IMatrix=new int[eqmat.NRows+1];
    for (i=0;i<eqmat.NRows+1;i++) {
      (*this).IMatrix[i]=eqmat.IMatrix[i];
    }
    if (eqmat.NonZeros>0) {
      (*this).Matrix=new double[eqmat.NonZeros];
      (*this).JMatrix=new int[eqmat.NonZeros];
      for (i=0;i<eqmat.NonZeros;i++) {
	(*this).Matrix[i]=eqmat.Matrix[i];
	(*this).JMatrix[i]=eqmat.JMatrix[i];
      }
    }
  }

  return *this;
}

CDReducedMatrix2 CDReducedMatrix2::operator+(CDReducedMatrix2 addmat)
{
  // This algorithm was taken from the book:
  // Sparse matrix technology , Author: Sergio Pissanetsky

  CDReducedMatrix2 res;
  vector<int> resimatrix,resjmatrix;
  vector<double> resmatrix;
  int *scratch,*resiptr,*resjptr;
  double *dscratch;
  int i,j,ip,jp,ih,iaa,iab,iba,ibb,ica,icb;

  if (NRows!=addmat.NRows) {
    printf("CDReducedMatrix2 CDReducedMatrix2::operator+ :\n");
    printf("Rows do not match! Aborting ..\n");
    exit(0);
  }
  if (NCols!=addmat.NCols) {
    printf("CDReducedMatrix2 CDReducedMatrix2::operator+ :\n");
    printf("Columns do not match! Aborting ..\n");
    exit(0);
  }

  res.ClearMatrix();
  resimatrix.clear();
  resjmatrix.clear();

  if ((NRows>0) && (NCols>0)) {
    // symbolical addition
    scratch=new int[NCols];
    ip=0;
    for (i=0;i<NCols;i++) { scratch[i]=-1; };
    for (i=0;i<NRows;i++) {
      resimatrix.push_back(ip);
      iaa=IMatrix[i];
      iab=IMatrix[i+1]-1;
      for (jp=iaa;jp<=iab;jp++) {
	j=JMatrix[jp];
	resjmatrix.push_back(j);
	ip=ip+1;
	scratch[j]=i;
      }
      iba=addmat.IMatrix[i];
      ibb=addmat.IMatrix[i+1]-1;
      for (jp=iba;jp<=ibb;jp++) {
	j=addmat.JMatrix[jp];
	if (scratch[j]!=i) {
	  resjmatrix.push_back(j);
	  ip=ip+1;
	}
      }
    }
    resimatrix.push_back(ip);
    delete [] scratch;

    // copy the nonzero structure to the result matrix
    res.NRows=NRows;
    res.NCols=NCols;
    res.NonZeros=resjmatrix.size();
    if (resimatrix.size()>0) {
      resiptr=&resimatrix.front();
      res.IMatrix=new int[NRows+1];
      for (i=0;i<NRows+1;i++) { res.IMatrix[i]=resiptr[i]; };
    }
    if (res.NonZeros>0) {
      resjptr=&resjmatrix.front();
      res.JMatrix=new int[res.NonZeros];
      res.Matrix=new double[res.NonZeros];
      for (i=0;i<res.NonZeros;i++) { res.JMatrix[i]=resjptr[i]; };
    }

    //    for (i=0;i<NRows+1;i++) { printf("%d ",res.IMatrix[i]); };
    //    printf("\n");

    // numerical addition
    dscratch=new double[NCols];
    for (i=0;i<NRows;i++) {
      ih=i+1;
      ica=res.IMatrix[i];
      icb=res.IMatrix[ih]-1;
      for (ip=ica;ip<=icb;ip++) { dscratch[res.JMatrix[ip]]=0; };
      iaa=IMatrix[i];
      iab=IMatrix[ih]-1;
      for (ip=iaa;ip<=iab;ip++) { dscratch[JMatrix[ip]]=Matrix[ip]; };
      iba=addmat.IMatrix[i];
      ibb=addmat.IMatrix[ih]-1;
      for (ip=iba;ip<=ibb;ip++) {
	j=addmat.JMatrix[ip];
	dscratch[j]=dscratch[j]+addmat.Matrix[ip];
      }
      for (ip=ica;ip<=icb;ip++) { res.Matrix[ip]=dscratch[res.JMatrix[ip]]; };
    }
    delete [] dscratch;
  }

  //  printf("This is not implemented yet! Aborting ..\n");
  //  exit(0);

  return res;
}

CDReducedMatrix2 CDReducedMatrix2::operator*(double factor)
{
  CDReducedMatrix2 res;
  int i;
  res=(*this);
  for (i=0;i<NonZeros;i++) { res.Matrix[i]=res.Matrix[i]*factor; };
  return res;
}

void CDReducedMatrix2::Mul(double factor)
{
  int i;
  for (i=0;i<NonZeros;i++) { Matrix[i]=Matrix[i]*factor; };
}

void CDReducedMatrix2::MulVec(double *invec,double *outvec)
{
  // Matrix vector multiplication.
  // 'invec' has to be an rray of size 'NCols' and must contain the input vector.
  // 'outvec' has to be an array of size 'NRows' and the will receive
  // the result vector.

  int i,j;

  for (i=0;i<NRows;i++) { outvec[i]=0.0; };
  for (i=0;i<NRows;i++) {
    for (j=IMatrix[i];j<IMatrix[i+1];j++) {
      outvec[i]=outvec[i]+Matrix[j]*invec[JMatrix[j]];
    }
  }
}

void CDReducedMatrix2::Mul3jAndPhase(int j1,int m1,int j2,int m2,
				     int j3,int m3,int phase)
{
  // Multiply the matrix with a wigner3j symbol and a phase.

  double factor;
  factor=wigner3j(j1,m1,j2,m2,j3,m3)*phase;
  Mul(factor);
}

void CDReducedMatrix2::Transpose(CDReducedMatrix2* transmat)
{
  // Transpose *this and give the result to 'transmat'.
  // *this will NOT be changed.
  // This routine runs in O(NonZeros).
  // Reference :Sparse matrix technology ; Author : Sergio Pissanetsky

  int i,j,icol;
  int *elecount;

  if (NonZeros==0) {
    transmat->ClearMatrix();
    transmat->NonZeros=0;
    transmat->NRows=NCols;
    transmat->NCols=NRows;
    return;
  }

  // Allocate arrays for the transposed matrix and set some NonZeros etc.
  transmat->NonZeros=NonZeros;
  transmat->NRows=NCols;
  transmat->NCols=NRows;
  transmat->IMatrix=new int[transmat->NRows+1];
  transmat->JMatrix=new int[transmat->NonZeros];
  transmat->Matrix=new double[transmat->NonZeros];

  // Count the number of elements of each column
  elecount=new int[NCols];

  for (i=0;i<NCols;i++) { elecount[i]=0; };

  for (i=0;i<NonZeros;i++) {
    j=JMatrix[i];
    elecount[j]=elecount[j]+1;
  }

  // test output
  //  for (i=0;i<NCols;i++) { printf("%d ",elecount[i]); };
  //  printf("\n");

  // The array elecount contains now in elecount[i] the number of nonzero
  // elements in column i. With that information it is easy to construct
  // IMatrix of the transposed matrix.
  transmat->IMatrix[0]=0;
  for (i=1;i<=transmat->NRows;i++) {
    transmat->IMatrix[i]=transmat->IMatrix[i-1]+elecount[i-1];
  }

  // test output
  //  for (i=0;i<=transmat->NRows;i++) { printf("%d ",transmat->IMatrix[i]); };
  //  printf("\n");

  // This loop does the real work. We use transmat->IMatrix as an array
  // of row pointers to remember at what position in transmat->Matrix we have
  // to include a new element. Every time we include a new element, we
  // raise the pointer. At the end we restore the transmat->IMatrix by using
  // 'elecount' as we did before.
  for (i=0;i<NRows;i++) { // i is now row-index of the original matrix
    for (j=IMatrix[i];j<IMatrix[i+1];j++) { // j is index og the elem. of '*this'
      icol=JMatrix[j]; // icol is the column index of the original matrix
      transmat->JMatrix[transmat->IMatrix[icol]]=i; // row of *this is new col.
      transmat->Matrix[transmat->IMatrix[icol]]=Matrix[j]; // copy element
      transmat->IMatrix[icol]=transmat->IMatrix[icol]+1; // raise pointer
    }
  }

  // test output
  //  for (i=0;i<=transmat->NRows;i++) { printf("%d ",transmat->IMatrix[i]); };
  //  printf("\n");

  // Restore transmat->IMatrix
  transmat->IMatrix[0]=0;
  for (i=1;i<=transmat->NRows;i++) {
    transmat->IMatrix[i]=transmat->IMatrix[i-1]+elecount[i-1];
  }

  delete [] elecount;
}

double CDReducedMatrix2::GetRMatElement(int* state1,int* state2,int optype,int ptnum,
					int* ptypes,CDidParticleCFP* cfptable)
{
  int i,phase;
  int J1,J2,J,J1p,J2p,Jp,k,j1,j2;
  int N1,N2,nu1,nu2,delta1,delta2;
  double w6j,cfp,res;

  /*
  printf("< ");
  for (i=0;i<ptnum*5;i++) { printf("%d ",state1[i]); };
  printf("|| %d || ",optype);
  for (i=0;i<ptnum*5;i++) { printf("%d ",state2[i]); };
  printf("> = ");
  */

  res=1.0;
  k=ptypes[abs(optype)-1]-1;

  // If ptnum>abs(optype) we apply:
  // <a1,J1,a2,J2,J||T^k(1)||a1p,J1p,a2p,J2p,Jp>_T =
  // (-)^{J1+J2+Jp+k}*sqrt((2*J+1)(2*Jp+1))*w6j(J1,J,J2,Jp,J1p,k)*
  // <a1,J1||T^k||a1p,J1p>_T * delta(a2,a2p)*delta(J2,J2p)
  // We don't have to check the delta functions. This routine assumes
  // these selections rules to be OK.
  while (ptnum>abs(optype)) {
    J1=state1[(ptnum-1)*5-1];
    J2=state1[ptnum*5-2];
    J=state1[ptnum*5-1];
    J1p=state2[(ptnum-1)*5-1];
    J2p=state2[ptnum*5-2];
    Jp=state2[ptnum*5-1];
    if ((J1+J2+Jp+k)%4==0) { phase=1; } else { phase=-1; };
    w6j=wigner6j(J1,J,J2,Jp,J1p,k);
    res=res*phase*sqrt((double)((J+1)*(Jp+1)))*w6j;
    ptnum=ptnum-1;
  }

  // Now we have abs(optype)==ptnum.
  N1=state1[ptnum*5-5];
  nu1=state1[ptnum*5-4];
  delta1=state1[ptnum*5-3];
  j1=state1[ptnum*5-2];
  N2=state2[ptnum*5-5];
  nu2=state2[ptnum*5-4];
  delta2=state2[ptnum*5-3];
  j2=state2[ptnum*5-2];
  if (optype>0) { // creator
    cfp=cfptable->GetCFP(N1,nu1,delta1,j1,nu2,delta2,j2);
    //      printf("+ %d %d %d %d %d %d %d %f %d\n",
    //	     N1,nu1,delta1,j1,nu2,delta2,j2,cfp,k);
    // convert the CFP to the reduced matrix element c+
    cfp=cfp*sqrt((double)(N1));
    // go to Talmi-convention of the red.mat.
    if ((2*j1+2*j2)%4==0) { phase=1; } else { phase=-1; };
    cfp=cfp*phase*sqrt((double)(j1+1));
    res=res*cfp;
  } else { // annihilator
    cfp=cfptable->GetCFP(N2,nu2,delta2,j2,nu1,delta1,j1);
    //      printf("~ %d %d %d %d %d %d %d %f %d\n",
    //	     N2,nu2,delta2,j2,nu1,delta1,j1,cfp,k);
    // convert the CFP to the reduced matrix element c+
    cfp=cfp*sqrt((double)(N2));
    // go to Talmi-convention of the red.mat.
    if ((2*j1+2*j2)%4==0) { phase=1; } else { phase=-1; };
    cfp=cfp*phase*sqrt((double)(j2+1));
    // convert the red. mat. element of c+ to the red. mat. of c~
    if ((j1-j2+k)%4==0) { phase=1; } else { phase=-1; };
    cfp=cfp*phase;
    res=res*cfp;
  }
  if (ptnum>1) {
    // apply:
    // <a1,J1,a2,J2,J||T^k(2)||a1p,J1p,a2p,J2p,Jp>_T =
    // (-)^{J1+J2p+J+k}*sqrt((2*J+1)(2*Jp+1))*w6j(J2,J,J1,Jp,J2p,k)*
    // <a2,J2||T^k||a2p,J2p>_T * delta(a1,a1p)*delta(J1,J1p)
    // The reduced matrix element of T^k was is already done.
    // This routine assumes that the delta functons are OK.
    // So there is left a phase, a square root and a wigner6j-symbol.
    J1=state1[(ptnum-1)*5-1];
    J2=state1[ptnum*5-2];
    J=state1[ptnum*5-1];
    J1p=state2[(ptnum-1)*5-1];
    J2p=state2[ptnum*5-2];
    Jp=state2[ptnum*5-1];
    if ((J1+J2p+J+k)%4==0) { phase=1; } else { phase=-1; };
    w6j=wigner6j(J2,J,J1,Jp,J2p,k);
    res=res*phase*sqrt((double)((J+1)*(Jp+1)))*w6j;
  }

  //  printf("%f\n",res);

  return res;
}

void CDReducedMatrix2::CalcReducedMatrix(CLinearBasisQN2* bas1,CLinearBasisQN2* bas2,
					 CDidParticleCFP* cfptable,int optype)
{
  int i,j,ic,indexcount;
  vector<int> indexlist,ptypes;
  vector<vector<int> > allilists;
  double redmat;

  ClearMatrix();
  ptypes=bas1->ptypes;

  // Get for each state of basis1 an indexlist of state of basis2.
  allilists.clear();
  indexcount=0;
  for (i=0;i<bas1->GetStateNum();i++) {
    bas2->SelectStates(bas1->GetPtr(i),optype,&indexlist);
    allilists.push_back(indexlist);
    indexcount=indexcount+indexlist.size();
    //    for (j=0;j<allilists[i].size();j++) { printf("%d ",allilists[i][j]); };
    //    printf("\n");
  }

  // Construct the sparse matrix with by using 'allilists'.
  NonZeros=indexcount;
  NRows=bas1->GetStateNum();
  NCols=bas2->GetStateNum();
  if (NonZeros>0) {
    Matrix=new double[NonZeros];
    JMatrix=new int[NonZeros];
    IMatrix=new int[NRows+1];

    IMatrix[0]=0;
    for (i=0;i<NRows;i++) {
      IMatrix[i+1]=IMatrix[i]+allilists[i].size();
    }

    ic=0;
    for (i=0;i<allilists.size();i++) {
      for (j=0;j<allilists[i].size();j++) {
	redmat=GetRMatElement(bas1->GetPtr(i),bas2->GetPtr(allilists[i][j]),
			      optype,ptypes.size(),&ptypes.front(),cfptable);
	JMatrix[ic]=allilists[i][j];
	Matrix[ic]=redmat;
	ic=ic+1;
      }
    }
  }
}

void CDReducedMatrix2::
CalcReducedMatrix(CDReducedMatrix2* rmat1,CDReducedMatrix2* rmat2,
		  int* js1,int* jsp,int* js2,int oj1,int oj2,int oj)
{
  // <a1,J1||[T^k1 x T^k2]^k||a2,J2>_T =
  // (-)^(J1+J2+k) * sqrt(2*k+1) *
  // sum_{ap,Jp} wigner6j(k1,k2,k,J2,J1,Jp} *
  // <a1,J1||T^k1||ap,Jp>_T * <ap,Jp||T^k2||a2,J2>_T
  //
  // Arguments: '*rmat1' = <a1,J1||T^k1||ap,Jp>_T
  //            '*rmat2' = <ap,Jp||T^k2||a2,J2>_T
  //            'js1'    = array with spins of the basis |a1,J1>
  //            'jsp'    = array with spins of the basis |ap,Jp>
  //            'js2'    = array with spins of the basis |a2,J2>
  //  'oj1','oj2','oj' are the spins of the operators in units of hbar/2

  int i,j,k;
  int J1,J2,Jp,phase;
  int indexsize,istart,length,jj,ibuf;
  int *index,*tempptr;
  double *temp;
  vector<int> tempjmatrix;
  double ajj,w6j;

  // Check the sizes of the matrices
  if (!(rmat1->NCols==rmat2->NRows)) {
    printf("void CDReducedMatrix2::CalcReducedMatrix(..) :\n");
    printf("Column number of the first matrix is different to row\n");
    printf("number of the second matrix! Aborting ..\n");
    exit(0);
  }

  if ((rmat1->NRows==0) || (rmat2->NCols==0) ||
      (rmat2->NRows==0) || (rmat2->NCols==0)) {
    SetNullMatrix(rmat1->NRows,rmat2->NCols);
    return;
  }

  if ((rmat1->NonZeros==0) || (rmat2->NonZeros==0)) {
    SetNullMatrix(rmat1->NRows,rmat2->NCols);
    return;
  }

  // This algorithm is basically a matrix multiplication. It works in two
  // steps. The first step creates the array IMatrix and JMatrix or in other
  // words, the first step determines the non zero structure of the result
  // matrix. The second step calculates the non zero matrix elements.
  // Reference :Sparse matrix technology ; Author : Sergio Pissanetsky

  ClearMatrix();

  tempjmatrix.clear();

  NRows=rmat1->NRows;
  NCols=rmat2->NCols;
  IMatrix=new int[NRows+1];

  // The array 'index' is used as a scratch. Lets initialize it.
  indexsize=rmat1->NRows;
  if (rmat2->NRows>indexsize) { indexsize=rmat2->NRows; };
  if (rmat2->NCols>indexsize) { indexsize=rmat2->NCols; };

  indexsize=indexsize+1;
  index=new int[indexsize];
  for (i=0;i<indexsize;i++) { index[i]=0; };

  // Now begin with the first step.
  IMatrix[0]=0;
  for (i=0;i<NRows;i++) { // i is row index of rmat1 and result matrix
    istart=-1;
    length=0;

    for (jj=rmat1->IMatrix[i];jj<rmat1->IMatrix[i+1];jj++) {
      j=rmat1->JMatrix[jj];

      for (k=rmat2->IMatrix[j];k<rmat2->IMatrix[j+1];k++) {
	//       	if (index[rmat2->JMatrix[k]]==0) { // bug ?
       	if (index[rmat2->JMatrix[k]+1]==0) {
	  //	  index[rmat2->JMatrix[k]]=istart; // bug ?
	  index[rmat2->JMatrix[k]+1]=istart;

	  //	  istart=rmat2->JMatrix[k]; // bug ?
	  istart=rmat2->JMatrix[k]+1;

	  length=length+1;
	}
      }
    }
    IMatrix[i+1]=IMatrix[i]+length;

    //    printf("debug : %d\n",length);

    for (j=IMatrix[i];j<IMatrix[i+1];j++) {
      //      printf("%d %d\n",j,istart);
      ibuf=istart;
      //      JMatrix[j]=istart;

      // tempjmatrix.push_back(ibuf); // bug ?
      if (ibuf>0) {
	tempjmatrix.push_back(ibuf-1);
      } else {
	tempjmatrix.push_back(ibuf);
      }

      istart=index[istart];
      //      index[JMatrix[j]]=0;
      index[ibuf]=0;
    }
    index[i]=0;
  }

  NonZeros=tempjmatrix.size();

  //  printf("debug : %d %d %d\n",NCols,NRows,NonZeros);

  if (NonZeros>0) {
    JMatrix=new int[NonZeros];
    Matrix=new double[NonZeros];
  }

  tempptr=&tempjmatrix.front();
  for (i=0;i<NonZeros;i++) { JMatrix[i]=tempptr[i]; };

  // First step is finished. Do the second step.

  // The array temp is used to store partial sums.
  temp=new double[indexsize];
  for (i=0;i<indexsize;i++) { temp[i]=0; };

  for (i=0;i<NRows;i++) {
    for (jj=rmat1->IMatrix[i];jj<rmat1->IMatrix[i+1];jj++) {
      j=rmat1->JMatrix[jj];
      ajj=rmat1->Matrix[jj];
      for (k=rmat2->IMatrix[j];k<rmat2->IMatrix[j+1];k++) {
	J1=js1[i];
	J2=js2[rmat2->JMatrix[k]];
	Jp=jsp[j];
	w6j=wigner6j(oj1,oj2,oj,J2,J1,Jp);
	temp[rmat2->JMatrix[k]]=temp[rmat2->JMatrix[k]]+w6j*ajj*rmat2->Matrix[k];
      }
    }
    for (j=IMatrix[i];j<IMatrix[i+1];j++) {
      Matrix[j]=temp[JMatrix[j]];
      temp[JMatrix[j]]=0;
    }
  }
  // Second step is finished.

  // Multiply each matrix element by a phase and a square root.
  for (i=0;i<NRows;i++) {
    for (j=IMatrix[i];j<IMatrix[i+1];j++) {
      k=JMatrix[j];
      if ((js1[i]+js2[k]+oj)%4==0) { phase=1; } else { phase=-1; };
      Matrix[j]=Matrix[j]*phase*sqrt((double)(oj+1));
    }
  }

  delete [] index;
  delete [] temp;
}

void CDReducedMatrix2::ClearMatrix()
{
  delete [] IMatrix;
  IMatrix=new int[1];
  IMatrix[0]=0;
  if (NonZeros>0) {
    delete [] Matrix;
    delete [] JMatrix;
  }
  NonZeros=0;
  NRows=0;
  NCols=0;
}

void CDReducedMatrix2::SetNullMatrix(int newrows,int newcols)
{
  int i;
  ClearMatrix();
  NonZeros=0;
  NRows=newrows;
  NCols=newcols;
  if (newrows+1>0) {
    IMatrix=new int[newrows+1];
    for (i=0;i<newrows+1;i++) { IMatrix[i]=0; };
  }
}

void CDReducedMatrix2::Print()
{
  int i,j,ic;
  double *row;

  printf("Number of non zero elements : %d\n",NonZeros);
  printf("Number of rows              : %d\n",NRows);
  printf("Number of columns           : %d\n",NCols);

  /*
  for (i=0;i<NRows+1;i++) { printf("%d ",IMatrix[i]); };
  printf("\n");
  for (i=0;i<NonZeros;i++) { printf("%d ",JMatrix[i]); };
  printf("\n");
  for (i=0;i<NonZeros;i++) { printf("%f ",Matrix[i]); };
  printf("\n");
  */

  if (NCols>0) {
    row=new double[NCols];

    for (i=0;i<NRows;i++) {
      for (j=0;j<NCols;j++) { row[j]=0; };
      for (j=IMatrix[i];j<IMatrix[i+1];j++) { row[JMatrix[j]]=Matrix[j]; };
      for (j=0;j<NCols;j++) { printf("% f ",row[j]); };
      printf("\n");
    }

    delete [] row;
  }
}

int CDReducedMatrix2::GetRowNum()
{
  return NRows;
}

int CDReducedMatrix2::GetColumnNum()
{
  return NCols;
}

void CDReducedMatrix2::GetMatrix(double *mat)
{
  // Copy the sparse matrix *this to an array of doubles. This array
  // must be already allocated and have the size = NCols*NRows.
  // The array 'mat' is NOT in any sparse matrix format. It will
  // get the first row, then the second and so on including all zeros.

  int i,j,k,i1,i2;

  // Clear the matrix
  k=0;
  for (i=0;i<NRows;i++) {
    for (j=0;j<NCols;j++) {
      mat[k]=0;
      k=k+1;
    }
  }

  // Now copy the nonzeros
  for (i=0;i<NRows;i++) {
    i1=IMatrix[i];
    i2=IMatrix[i+1];
    for (j=i1;j<i2;j++) {
      k=JMatrix[j];
      mat[i*NCols+k]=Matrix[j];
    }
  }
}

void CDReducedMatrix2::WriteToFile(char *filename)
{
  int i;
  ofstream fout(filename,ios::out);

  if (!fout.good()) {
    printf("void CDReducedMatrix2::WriteToFile(char *filename) :\n");
    printf("Cannot open file for writing. Aborting ..\n");
    exit(0);
  }

  fout.precision(20); // should be just enough

  fout << NonZeros << endl;
  fout << NRows << endl;
  fout << NCols << endl;
  for (i=0;i<NRows+1;i++) { fout << IMatrix[i] << " "; };
  fout << endl;
  for (i=0;i<NonZeros;i++) { fout << JMatrix[i] << " "; };
  fout << endl;
  for (i=0;i<NonZeros;i++) { fout << Matrix[i] << " "; };
  fout << endl;

  fout.close();
}

////////////////////////////////////////////////////////////////
//////////////////// PRIVATE FUNCTIONS /////////////////////////
////////////////////////////////////////////////////////////////

double CDReducedMatrix2::wigner6j(int J1,int J2,int J3,int L1,int L2,int L3)
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

double CDReducedMatrix2::wigner6jhelp(int J1,int J2,int J3)
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
    cout << "double wigner6jhelp(int,int,int) :" << endl;
    cout << "Index out of array range, you have to ";
    cout << "enlarge FCT[] !" << endl;
    exit(0);
  }
  FDELTA=sqrt(FCT[IW1]*FCT[IW2]*FCT[IW3]/FCT[IW4]);
  return FDELTA/sqrt((double)10.0);
}

double CDReducedMatrix2::clebsch(int J1,int M1,int J2,int M2,int J,int M)
{
  // all in units hbar/2 !!!
  // gives back <J1,M1,J2,M2|J,M>

  long Z,ZMIN,ZMAX,FASE,F1,F2;
  long JA,JB,JC,JD,JE,JF,JG,JH,JI,JJ;
  double CC;
  CC=0;
  if (M1+M2==M) {
    if (abs(M1)<=abs(J1)) {
      if (abs(M2)<=abs(J2)) {
        if (abs(M)<=abs(J)) {
          if (J<=J1+J2) {
            if (J>=abs(J1-J2)) {
              ZMIN=0;
              if (J-J2+M1<0) {
                ZMIN=-J+J2-M1;
              }
              if (J-J1-M2+ZMIN<0) {
                ZMIN=-J+J1+M2;
              }
              ZMAX=J1+J2-J;
              if (J2+M2-ZMAX<0) {
                ZMAX=J2+M2;
              }
              if (J1-M1-ZMAX<0) {
                ZMAX=J1-M1;
              }
              for (Z=ZMIN;Z<=ZMAX;Z=Z+2) {
                JA=Z/2+1;
                JB=(J1+J2-J-Z)/2+1;
                JC=(J1-M1-Z)/2+1;
                JD=(J2+M2-Z)/2+1;
                JE=(J-J2+M1+Z)/2+1;
                JF=(J-J1-M2+Z)/2+1;
                FASE=(long)pow((double)(-1),(double)(Z/2));
                F2=FASE;
		if ((JA>=FCT_Size) || (JB>=FCT_Size) || (JC>=FCT_Size) ||
		    (JD>=FCT_Size) || (JE>=FCT_Size) || (JF>=FCT_Size)) {
		  // index out of array size , exit
		  cout << "double clebsch";
		  cout << "(int,int,int,int,int,int) :" << endl;
		  cout << "Index out of array range, you have to ";
		  cout << "enlarge FCT[] !" << endl;
		  exit(0);
		}
                CC=CC+F2/(FCT[JA]*FCT[JB]*FCT[JC]*FCT[JD]*FCT[JE]*FCT[JF]);
              }
              JA=(J1+J2-J)/2+1;
              JB=(J1-J2+J)/2+1;
              JC=(-J1+J2+J)/2+1;
              JD=(J1+M1)/2+1;
              JE=(J1-M1)/2+1;
              JF=(J2+M2)/2+1;
              JG=(J2-M2)/2+1;
              JH=(J+M)/2+1;
              JI=(J-M)/2+1;
              JJ=(J1+J2+J+2)/2+1;
              F1=J+1;
	      if ((JA>=FCT_Size) || (JB>=FCT_Size) || (JC>=FCT_Size) ||
		  (JD>=FCT_Size) || (JE>=FCT_Size) || (JF>=FCT_Size) ||
		  (JG>=FCT_Size) || (JH>=FCT_Size) || (JI>=FCT_Size) ||
		  (JJ>=FCT_Size)) {
		// index out of array size , exit
		cout << "double clebsch";
		cout << "(int,int,int,int,int,int) :" << endl;
		cout << "Index out of array range, you have to ";
		cout << "enlarge FCT[] !" << endl;
		exit(0);
	      }
              CC=sqrt((double)(F1*FCT[JA]*FCT[JB]*FCT[JC]*FCT[JD]
			      *FCT[JE]*FCT[JF]*FCT[JG]*FCT[JH]*
			      FCT[JI]/FCT[JJ]))*CC;
            }
          }
        }
      }
    }
  }
  CC=CC/sqrt((double)10);
  if (CC*CC<10e-20) {
    CC=0;
  }
  return CC;
}

double CDReducedMatrix2::wigner3j(int j1,int m1,int j2,int m2,int j3,int m3)
{
  // all in units of hbar/2
  // / j1 j2 J \ = \frac{(-1)^{j1-j2-M}}{\sqrt{2*J+1}} * <j1,m1,j2,m2|J,-M>
  // \ m1 m2 M /
  double res;

  res=clebsch(j1,m1,j2,m2,j3,-m3)*pow((double)(-1),
				      (double)((j1-j2-m3)/2.0))/
    sqrt((double)(j3+1));

  return res;
}
