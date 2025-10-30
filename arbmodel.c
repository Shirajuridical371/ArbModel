///////////////////////////////////////////////////

#include <cstdlib>
#include <string>
#include <fstream>
#include <iostream>

#include "./Classes/clinearbasisqn2.h"
#include "./Classes/cdreducedmatrix2.h"
#include "./Classes/cdidparticlecfp.h"
#include "./Classes/cdrecouple.h"
#include "./Classes/csubroutines.h"
#include "./Classes/carpackwrapper.h"

////////////////////////////////////////////////////
int model_ptnum;
string model_isf_dir;
vector<string> model_ptnames_array;
map<string,int> model_ptnames;
vector<int> model_ptparities;
vector<int> model_ptsymmetries;
vector<int> model_ptmaxnums;
vector<int> model_basislimits;
vector<int> model_filterseniorities;
int model_totalparticlenum;
vector<int> model_totalspin;
vector<int> model_totalparity;
vector<int> model_states_to_calc;
int model_groundstateindex;
vector<int> model_usehamiltonian;
vector<int> model_linehamilt;
vector<string> model_printbasis;
vector<string> model_printeigenvals;
vector<string> model_printeigenvecs;
string model_printinputsummary;
int model_verboselevel;
vector<int> model_diag_method;
vector<string> model_hamilttofile;

vector<string> hamilt_names_array;
map<string,int> hamilt_names;
vector<CDRecouple> hamiltonian;

// For each spin and parity we have one hamiltonian matrix, one set of
// eigenvalues and a set of eigenvectors and a basis in terms of quantum numbers
vector<CDReducedMatrix2> allhsparsemats; // to store all the hamiltonian matrices
vector<int> allhsizes;                  // to store all hamiltonian sizes
vector<vector<double> > alleigenvals;   // to store all eigenvalues
vector<vector<double> > alleigenvecs;   // to store all eigenvectors
vector<CLinearBasisQN2> allbasis;        // to store all basis quantum numbers

// For storing some one body operators (for example TE2)
vector<string> allobopsnames_array;    // all names of the operators
map<string,int> allobopsnames;         // for mapping names <-> indexes
vector<vector<int> > allobops;         // the operators and the
vector<vector<double> > allobopscoefs; // corresponding coefficients

// Array for storing which matrix elements should be calculated.
// Five integers for each matrix element.
// 1. Reference number of the calculation for the bra state (starts with 1)
// 2. Number of the bra state (starts wit 1)
// 3. index of the operator name (for array 'allobopsnames_array, starts with 0)
// 4. Reference number of the calculation for the ket state (starts with 1)
// 5. Number of the ket state (starts wit 1)
vector<int> matels_to_calc;

// For each matrix element to calculate we have one pointer to an operator
// matrix. By using this array of pointers, we avoid storing a matrix double.
vector<int> allobsparsematsptr; // to store the indices (pointers)
vector<CDReducedMatrix2> allobsparsemats; // to store all one-body operator matrices
vector<double> allobmatels;
////////////////////////////////////////////////////

void clear_parameters()
{
  model_isf_dir=string("");
  model_ptnum=0;
  model_ptnames_array.clear();
  model_ptnames.clear();
  model_ptparities.clear();
  model_ptsymmetries.clear();
  model_ptmaxnums.clear();
  model_basislimits.clear();
  model_filterseniorities.clear();
  model_totalparticlenum=0;
  model_totalspin.clear();
  model_totalparity.clear();
  model_states_to_calc.clear();
  model_groundstateindex=-1;
  model_usehamiltonian.clear();
  model_linehamilt.clear();
  model_printbasis.clear();
  model_printeigenvals.clear();
  model_printeigenvecs.clear();
  model_printinputsummary=string("y");
  model_verboselevel=1;
  model_diag_method.clear();
  model_hamilttofile.clear();

  hamilt_names_array.clear();
  hamilt_names.clear();
  hamiltonian.clear();

  allobopsnames_array.clear();
  allobopsnames.clear();
  allobops.clear();
  allobopscoefs.clear();
  matels_to_calc.clear();
}

void print_input_summary()
{
  int i,count;
  map<string,int>::iterator nameiter;

  printf("Folder of the isoscalara factors : %s\n",model_isf_dir.c_str());

  printf("Number of distinguishable particles : %d\n",model_ptnum);
  printf("Names of the particles              : ");
  for (i=0;i<model_ptnames_array.size();i++) {
    printf("%s ",model_ptnames_array[i].c_str());
  }
  printf("\n");
  printf("Parities of the particles           : ");
  for (i=0;i<model_ptparities.size();i++) { printf("% d ",model_ptparities[i]); };
  printf("\n");
  printf("Symmetriesies of the particles      : ");
  for (i=0;i<model_ptsymmetries.size();i++) {
    printf("U(%d) ",model_ptsymmetries[i]);
  }
  printf("\n");
  printf("Maximum number of the particles     : ");
  for (i=0;i<model_ptmaxnums.size();i++) { printf("%d ",model_ptmaxnums[i]); };
  printf("\n");
  printf("Seniorities of the basis states     : ");
  for (i=0;i<model_filterseniorities.size();i++) {
    if (model_filterseniorities[i]==-1) {
      printf("* ");
    } else {
      printf("%d ",model_filterseniorities[i]);
    }
  }
  printf("\n");
  for (i=0;i<model_basislimits.size();i++) {
    if (i%(model_ptnum+2)==0) { printf("Basis limits : "); };
    if (i%(model_ptnum+2)==model_ptnum) {
      if (model_basislimits[i]==1) { printf("= "); };
      if (model_basislimits[i]==2) { printf("<= "); };
    } else {
      if (i%(model_ptnum+2)==model_ptnum+1) {
	printf("%d ",model_basislimits[i]);
      } else {
	if (i%(model_ptnum+2)==model_ptnum-1) {
	  printf("%d*n_%d ",model_basislimits[i],i%(model_ptnum+2)+1);
	} else {
	  printf("%d*n_%d + ",model_basislimits[i],i%(model_ptnum+2)+1);
	}
      }
    }
    if (i%(model_ptnum+2)==model_ptnum+1) { printf("\n"); };
  }
  printf("Total particle number of the states : %d\n",model_totalparticlenum);
  printf("Reference number of the calculation : ");
  for (i=0;i<model_totalspin.size();i++) { printf("%3d ",i+1); };
  printf("\n");
  printf("Total spin(s) of the states         : ");
  for (i=0;i<model_totalspin.size();i++) { printf("%3d ",model_totalspin[i]); };
  printf("\n");
  printf("Total parity(s) of the states       : ");
  for (i=0;i<model_totalparity.size();i++) { printf("%3d ",model_totalparity[i]); };
  printf("\n");
  if (model_groundstateindex==-1) {
    printf("No groundstate was given.\n");
  } else {
    printf("Groundstate is the first state of the %d-th calculation\n",
	   model_groundstateindex+1);
  }
  printf("Number of states to calculate       : ");
  for (i=0;i<model_states_to_calc.size();i++) {
    printf("%3d ",model_states_to_calc[i]);
  }
  printf("\n");
  printf("Use Hamiltonian                     : ");
  for (i=0;i<model_usehamiltonian.size();i++) {
    printf("%s ",hamilt_names_array[model_usehamiltonian[i]].c_str());
  }
  printf("\n");
  printf("Diagonalisation method              : ");
  for (i=0;i<model_diag_method.size();i++) { printf("%3d ",model_diag_method[i]); };
  printf("\n");
  printf("Files for storing the hamiltonian   : ");
  for (i=0;i<model_hamilttofile.size();i++) {
    printf("%s ",model_hamilttofile[i].c_str());
  }
  printf("\n");  

  if (model_linehamilt.size()==3) {
    printf("Line Hamiltonian                    : %d %d %d\n",model_linehamilt[0],
	   model_linehamilt[1],model_linehamilt[2]);
  } else {
    printf("Doing no line hamiltonian.\n");
  }
  printf("Print out the basis vectors         : ");
  for (i=0;i<model_printbasis.size();i++) {
    printf("  %s ",model_printbasis[i].c_str());
  }
  printf("\n");
  printf("Print out the eigenvalues           : ");
  for (i=0;i<model_printeigenvals.size();i++) {
    printf("  %s ",model_printeigenvals[i].c_str());
  }
  printf("\n");
  printf("Print out the eigenvectors          : ");
  for (i=0;i<model_printeigenvecs.size();i++) {
    printf("  %s ",model_printeigenvecs[i].c_str());
  }
  printf("\n");
  printf("Print input summary                 : %s\n",model_printinputsummary.c_str());
  printf("Verbose level                       : %d\n",model_verboselevel);

  printf("Number of defined operators : %d\n",(int)allobopsnames.size());
  for (count=0;count<allobopsnames_array.size();count++) {
    printf("%s = ",allobopsnames_array[count].c_str());
    for (i=0;i<allobops[count].size()/3;i++) {
      if (i!=0) { printf(" + "); };
      printf("%e * [%s+ x %s~]^(%d/2)",allobopscoefs[count][i],
	     model_ptnames_array[allobops[count][3*i]-1].c_str(),
	     model_ptnames_array[allobops[count][3*i+1]-1].c_str(),
	     allobops[count][3*i+2]);
    }
    printf("\n");
  }

  printf("Matrix elements to calculate :\n");
  for (i=0;i<matels_to_calc.size()/5;i++) {
    printf("<%d %d || %s || %d %d>\n",matels_to_calc[i*5],matels_to_calc[i*5+1],
	   allobopsnames_array[matels_to_calc[i*5+2]].c_str(),
	   matels_to_calc[i*5+3],matels_to_calc[i*5+4]);
  }

  printf("Using the hamiltonian(s) :\n");
  for (i=0;i<hamiltonian.size();i++) {
    printf("%s :\n",hamilt_names_array[i].c_str());
    hamiltonian[i].Print();
  }
}

void check_parameter()
{
  int i,j;
  if (model_isf_dir==string("")) {
    printf("No folder where to find the isoscalar factors was given! Aborting ..\n");
    exit(0);
  }
  if (!((model_ptnum>=1) && (model_ptnum<=10))) {
    printf("Number of distinguishable particles out of range! Aborting ..\n");
    exit(0);
  }
  if (model_ptnames.size()!=model_ptnum) {
    printf("Wrong number of particle names! Aborting..\n");
    exit(0);
  }
  if (model_ptparities.size()!=model_ptnum) {
    printf("Wrong number of particle parities! Aborting..\n");
    exit(0);
  }
  if (model_ptsymmetries.size()!=model_ptnum) {
    printf("Wrong number of particle symmetries! Aborting..\n");
    exit(0);
  }
  if (model_ptmaxnums.size()!=model_ptnum) {
    printf("Wrong number of maximum particle numbers! Aborting..\n");
    exit(0);
  }
  if (model_basislimits.size()%(model_ptnum+2)!=0) {
    printf("There is something wrong with the basis limits. Aborting ..\n");
    exit(0);
  }
  if (model_filterseniorities.size()!=model_ptnum) {
    printf("Wrong number of seniorities for the filter function. Aborting ..\n");
    exit(0);
  }
  for (i=0;i<model_filterseniorities.size();i++) {
    if (model_filterseniorities[i]<-1) {
      printf("Negative seniority given for the filter function. Aborting ..\n");
      exit(0);
    }
  }
  if (model_totalparticlenum<=1) {
    printf("Particle number has to be >= 2. Aborting ..\n");
    exit(0);
  }
  //  if (!((model_totalparticlenum>=3) && (model_totalparticlenum<=25))) {
  //    printf("Total particle number out of range! Aborting ..\n");
  //    exit(0);
  //  }
  for (i=0;i<model_totalspin.size();i++) {
    if (!((model_totalspin[i]>=0) && (model_totalspin[i]<=64))) {
      printf("Total spin out of range! Aborting ..\n");
      exit(0);
    }
  }
  for (i=0;i<model_totalparity.size();i++) {
    if (!((model_totalparity[i]==1) || (model_totalparity[i]==-1))) {
      printf("Total parity out of range! Aborting ..\n");
      exit(0);
    }
  }
  if (model_totalspin.size()!=model_totalparity.size()) {
    printf("Number of total spins does not match the number of\n");
    printf("total parities! Aborting ..\n");
    exit(0);
  }
  if (model_totalspin.size()!=model_states_to_calc.size()) {
    printf("Numbers of the states to calculate does not match the number\n");
    printf("of total spins! Aborting ..\n");
    exit(0);
  }
  if (model_totalspin.size()!=model_printbasis.size()) {
    printf("Number of flags for printing out the basis does not match\n");
    printf("the number of total spins! Aborting ..\n");
    exit(0);
  }
  if (model_totalspin.size()!=model_printeigenvals.size()) {
    printf("Number of flags for printing out the eigenvalues does not match\n");
    printf("the number of total spins! Aborting ..\n");
    exit(0);
  }
  if (model_totalspin.size()!=model_printeigenvecs.size()) {
    printf("Number of flags for printing out the eigenvectors does not match\n");
    printf("the number of total spins! Aborting ..\n");
    exit(0);
  }
  if (model_totalspin.size()!=model_usehamiltonian.size()) {
    printf("Number of the hamiltonians to use does not match\n");
    printf("the number of total spins! Aborting ..\n");
    exit(0);
  }
  if (model_totalspin.size()!=model_diag_method.size()) {
    printf("Number of diagonalisation methods to use does not match\n");
    printf("the number of hamiltonians to use! Aborting ..\n");
    exit(0);
  }
  if ((model_hamilttofile.size()!=0) &&
      (model_hamilttofile.size()!=model_totalspin.size())) {
    printf("The number of filenames for writing the hamiltonian matrix to has\n");
    printf("to be either zero or muat match the number of hamiltonians! ");
    printf("Aborting ..\n");
    exit(0);
  }
  for (i=0;i<model_diag_method.size();i++) {
    j=model_diag_method[i];
    if (!( (j>=1) && (j<=1) ) ) {
      printf("Unknown diagonalisation method : %d ! Aborting ..\n",j);
      exit(0);
    }
  }

  for (i=0;i<matels_to_calc.size()/5;i++) {
    if (!((matels_to_calc[i*5+1]<=model_states_to_calc[matels_to_calc[i*5]-1]) ||
	  (model_states_to_calc[matels_to_calc[i*5]-1]==-1))) {
      printf("NumberOfStatesToCalc contradicts with the requested matrix elements!\n");
      printf("Aborting ..\n");
      exit(0);
    }
    if (!((matels_to_calc[i*5+4]<=model_states_to_calc[matels_to_calc[i*5+3]-1]) ||
	  (model_states_to_calc[matels_to_calc[i*5+3]-1]==-1))) {
      printf("NumberOfStatesToCalc contradicts with the requested matrix elements!\n");
      printf("Aborting ..\n");
      exit(0);
    }
  }
  if (!((model_verboselevel>=0) && (model_verboselevel<=2))) {
    printf("Verbose level out of range. Aborting ..\n");
    exit(0);
  }
  if (model_linehamilt.size()==3) { // if we do a line hamiltonian
    if (!((model_linehamilt[0]>=1) && (model_linehamilt[0]<=model_totalspin.size()) &&
	  (model_linehamilt[1]>=1) && (model_linehamilt[1]<=model_totalspin.size()))) {
      printf("Reference numbers out of range for doing a line hamiltonian!\n");
      printf("Aborting ..\n");
      exit(0);
    }
    if (model_totalspin[model_linehamilt[0]-1]!=
	model_totalspin[model_linehamilt[1]-1]) {
      printf("Total spins for the line hamiltonian are different! Aborting ..\n");
      exit(0);
    }
    if (model_totalparity[model_linehamilt[0]-1]!=
	model_totalparity[model_linehamilt[1]-1]) {
      printf("Total parities for the line hamiltonian are different! Aborting ..\n");
      exit(0);
    }
    if (model_states_to_calc[model_linehamilt[0]-1]!=
	model_states_to_calc[model_linehamilt[1]-1]) {
      printf("Number of states to calculate are different for both sides\n");
      printf("of the line hamiltonian. Abrting ..\n");
      exit(0);
    }
    if (model_diag_method[model_linehamilt[0]-1]!=
	model_diag_method[model_linehamilt[1]-1]) {
      printf("Diagonalisation methods are different for both sides\n");
      printf("of the line hamiltonian. Aborting ..\n");
      exit(0);
    }
    if (!(model_linehamilt[2]>1)) {
      printf("Number of points for the line hamiltonian is < 2. Aborting...\n");
      exit(0);
    }
  }
  if (model_groundstateindex!=-1) {
    if (model_groundstateindex!=0) {
      printf("A groundstate was given, but the number of the corresponding\n");
      printf("calculation is out of range!\n");
      printf("Aborting ..\n");
      exit(0);
    }
  }
}

void process_hamilt_line(vector<string> linedata,double factor,int* lineok)
{
  map<string,int>::iterator nameiter;
  char delimiter='\n';
  string inputline,inputstr;
  vector<string> linedata2;
  int h1,h2,h3,h4,h5;
  double hcoef;
  int flag,linecount,lineok2,hindex,i;
  CDRecouple emptyhamilt;

  emptyhamilt.ClearAll();

  // Check if we already have a hamiltonian with the given name.
  if (linedata.size()>=2) {
    nameiter=hamilt_names.find(linedata[1]);
    if (nameiter==hamilt_names.end()) { // new hamiltonian
      hamilt_names[linedata[1]]=hamiltonian.size();
      hamilt_names_array.push_back(linedata[1]);
      hamiltonian.push_back(emptyhamilt);
      (hamiltonian.end()-1)->SetParticleTypes(model_ptsymmetries);
      hindex=hamiltonian.size()-1;
    } else {
      hindex=nameiter->second;
    }
  }
  // Now add th term
  if (linedata.size()==7) { // one body term
    flag=1;
    nameiter=model_ptnames.find(linedata[2]);
    if (nameiter!=model_ptnames.end()) {
      h1=nameiter->second;
    } else { flag=0; };
    nameiter=model_ptnames.find(linedata[4]);
    if (nameiter!=model_ptnames.end()) {
      h2=nameiter->second;
    } else { flag=0; };
    hcoef=atof(linedata[6].c_str());
    if (flag==1) {
      if ((linedata[3]==string("+")) && (linedata[5]==string("~"))) {
	hamiltonian[hindex].Add_onebody(h1,h2,hcoef*factor);
	*lineok=1;
      }
    }
  }
  if (linedata.size()==12) { // two body term
    flag=1;
    nameiter=model_ptnames.find(linedata[2]);
    if (nameiter!=model_ptnames.end()) {
      h1=nameiter->second;
    } else { flag=0; };
    nameiter=model_ptnames.find(linedata[4]);
    if (nameiter!=model_ptnames.end()) {
      h2=nameiter->second;
    } else { flag=0; };
    nameiter=model_ptnames.find(linedata[6]);
    if (nameiter!=model_ptnames.end()) {
      h3=nameiter->second;
    } else { flag=0; };
    nameiter=model_ptnames.find(linedata[8]);
    if (nameiter!=model_ptnames.end()) {
      h4=nameiter->second;
    } else { flag=0; };
    h5=atoi(linedata[10].c_str());
    hcoef=atof(linedata[11].c_str());
    if (flag==1) {
      if ((linedata[3]==string("+")) && (linedata[5]==string("+")) &&
	  (linedata[7]==string("~")) && (linedata[9]==string("~"))) {
	hamiltonian[hindex].Add_ddtt_to_ddtt(h1,h2,h3,h4,h5,hcoef*factor);
	*lineok=1;
      }
      if ((linedata[3]==string("+")) && (linedata[5]==string("~")) &&
	  (linedata[7]==string("+")) && (linedata[9]==string("~"))) {
	hamiltonian[hindex].Add_dtdt_to_ddtt(h1,h2,h3,h4,h5,hcoef*factor);
	*lineok=1;
      }
    }
  }
}

void read_inputfile(char* filename)
{
  string inputline,inputstr;
  std::ifstream fin(filename, std::ios::in);
  char delimiter='\n';
  vector<string> linedata;
  int i,flag,lineok,linecount;
  map<string,int>::iterator nameiter;

  vector<int> newop;
  vector<double> newopcoefs;

  int h1,h2,h3,h4,h5,p1,p2,angmom;
  double hcoef,coef;

  vector<string> tmphnames;

  tmphnames.clear();

  if (!fin.is_open()) {
    std::cerr << "Cannot open file: " << filename << ". Aborting ..\n";
    std::exit(EXIT_FAILURE);
  }

  clear_parameters();

  if (fin.good()) {
    getline(fin,inputline,delimiter);
  }
  linecount=1;
  while (fin.good()) {
    //    cout << inputline << endl;
    // We have a complete line in 'inputline'.
    // Now we split this line at the space characters and put each element
    // into 'linedata'.
    linedata.clear();
    inputstr.clear();
    flag=0;
    for (i=0;i<inputline.length();i++) {
      if (flag==0) {
	if (!isspace(inputline[i])) {
	  inputstr.push_back(inputline[i]);
	  flag=1;
	}
      } else {
	if (!isspace(inputline[i])) {
	  inputstr.push_back(inputline[i]);
	} else {
	  linedata.push_back(inputstr);
	  inputstr.clear();
	  flag=0;
	}
      }
    }
    if (inputstr.length()>0) { linedata.push_back(inputstr); };
    //    for (i=0;i<linedata.size();i++) { cout << linedata[i] << endl; };

    if (linedata.size()>0) { // do we have something new? (maybe the line was empty)
      if (linedata[0][0]=='#') {
	lineok=1;
      } else {
	lineok=0;

	if (linedata[0]==string("FolderOfISF")) {
	  if (linedata.size()==2) {
	    model_isf_dir=linedata[1];
	    lineok=1;
	  }
	}

	if (linedata[0]==string("ParticleTypeNumber")) {
	  if (linedata.size()==2) { 
	    model_ptnum=atoi(linedata[1].c_str());
	    lineok=1;
	  }
	}
	if (linedata[0]==string("ParticleTypeNames")) {
	  if (linedata.size()>1) {
	    for (i=1;i<linedata.size();i++) {
	      model_ptnames_array.push_back(linedata[i]);
	      model_ptnames[linedata[i]]=i;
	    }
	    lineok=1;
	  }
	}
	if (linedata[0]==string("ParticleTypeParities")) {
	  if (linedata.size()>1) {
	    for (i=1;i<linedata.size();i++) {
	      model_ptparities.push_back(atoi(linedata[i].c_str()));
	    }
	    lineok=1;
	  }
	}
	if (linedata[0]==string("ParticleTypeSymmetries")) {
	  if (linedata.size()>1) {
	    for (i=1;i<linedata.size();i++) {
	      model_ptsymmetries.push_back(atoi(linedata[i].c_str()));
	    }
	    lineok=1;
	  }
	  //	  hamiltonian.SetParticleTypes(model_ptsymmetries);
	  for (i=0;i<hamiltonian.size();i++) {
	    hamiltonian[i].SetParticleTypes(model_ptsymmetries);
	  }
	}
	if (linedata[0]==string("ParticleTypeMaxNums")) {
	  if (linedata.size()>1) {
	    for (i=1;i<linedata.size();i++) {
	      model_ptmaxnums.push_back(atoi(linedata[i].c_str()));
	    }
	    lineok=1;
	  }
	}
	if (linedata[0]==string("ParticleTypeSeniorities")) {
	  if (linedata.size()>1) {
	    for (i=1;i<linedata.size();i++) {
	      model_filterseniorities.push_back(atoi(linedata[i].c_str()));
	    }
	    lineok=1;
	  }
	}
	if (linedata[0]==string("BasisLimit")) {
	  if (model_ptnum==0) {
	    printf("You cannot give a 'BasisLimit' before 'ParticleTypeNumber'.\n");
	  } else {
	    if (linedata.size()==model_ptnum+3) {
	      for (i=1;i<linedata.size();i++) {
		if (linedata[i]==string("=")) {
		  model_basislimits.push_back(1);
		} else {
		  if (linedata[i]==string("<=")) {
		    model_basislimits.push_back(2);
		  } else {
		    model_basislimits.push_back(atoi(linedata[i].c_str()));
		  }
		}
	      }
	      lineok=1;
	    }
	  }
	}
	if (linedata[0]==string("TotalParticleNumber")) {
	  if (linedata.size()==2) { 
	    model_totalparticlenum=atoi(linedata[1].c_str());
	    lineok=1;
	  }
	}
	if (linedata[0]==string("TotalSpin")) {
	  if (linedata.size()>=2) { 
	    for (i=1;i<linedata.size();i++) {
	      model_totalspin.push_back(atoi(linedata[i].c_str()));
	    }
	    lineok=1;
	  }
	}
	if (linedata[0]==string("TotalParity")) {
	  if (linedata.size()>=2) {
	    for (i=1;i<linedata.size();i++) {
	      model_totalparity.push_back(atoi(linedata[i].c_str()));
	    }
	    lineok=1;
	  }
	}
	if (linedata[0]==string("Groundstate")) {
	  if (linedata.size()==2) {
	    model_groundstateindex=atoi(linedata[1].c_str())-1;
	    lineok=1;
	  }
	}
	if (linedata[0]==string("H")) {
	  process_hamilt_line(linedata,1.0,&lineok);
	}
	if (linedata[0]==string("O")) {
	  if (((linedata.size()-2)%6==0) && (linedata.size()>=8)) {
	    nameiter=allobopsnames.find(linedata[1]);
	    if (nameiter==allobopsnames.end()) { // new operator
	      allobopsnames[linedata[1]]=allobops.size();
	      allobopsnames_array.push_back(linedata[1]);
	      newop.clear();
	      newopcoefs.clear();
	      lineok=1;
	      for(i=0;i<(linedata.size()-2)/6;i++) {
		nameiter=model_ptnames.find(linedata[2+i*6]);
		if (nameiter==model_ptnames.end()) {
		  printf("Unknown creator!\n");
		  lineok=0;
		}
		p1=nameiter->second;
		nameiter=model_ptnames.find(linedata[2+i*6+2]);
		if (nameiter==model_ptnames.end()) {
		  printf("Unknown annihilator!\n");
		  lineok=0;
		}
		p2=nameiter->second;
		if (!((linedata[2+i*6+1]==string("+")) &&
		      (linedata[2+i*6+3]==string("~")))) {
		  lineok=0;
		}
		if (lineok==1) {
		  newop.push_back(p1);
		  newop.push_back(p2);
		  newop.push_back(atoi(linedata[2+i*6+4].c_str()));
		  newopcoefs.push_back(atof(linedata[2+i*6+5].c_str()));
		}
	      }
	      if (lineok==1) {
		allobops.push_back(newop);
		allobopscoefs.push_back(newopcoefs);
	      }
	    } else { // operator does already exist
	      printf("Redefinition of the operator '%s'.\n",linedata[1].c_str());
	    }
	  }
	}
	if (linedata[0]==string("M")) {
	  if (linedata.size()==6) {
	    nameiter=allobopsnames.find(linedata[3]);
	    if (nameiter!=allobopsnames.end()) {
	      matels_to_calc.push_back(atoi(linedata[1].c_str()));
	      matels_to_calc.push_back(atoi(linedata[2].c_str()));
	      matels_to_calc.push_back(nameiter->second);
	      matels_to_calc.push_back(atoi(linedata[4].c_str()));
	      matels_to_calc.push_back(atoi(linedata[5].c_str()));
	      lineok=1;
	    } else {
	      printf("Unknown operator.\n");
	    }
	  }
	}
	if (linedata[0]==string("NumberOfStatesToCalc")) {
	  if (linedata.size()>=2) { 
	    for (i=1;i<linedata.size();i++) {
	      model_states_to_calc.push_back(atoi(linedata[i].c_str()));
	    }
	    lineok=1;
	  }
	}
	if (linedata[0]==string("UseHamiltonian")) {
	  if (linedata.size()>=2) {
	    tmphnames.clear();
	    for (i=1;i<linedata.size();i++) {
	      tmphnames.push_back(linedata[i].c_str());
	    }
	    lineok=1;
	  }
	}
	if (linedata[0]==string("DiagMethod")) {
	  if (linedata.size()>=2) {
	    for (i=1;i<linedata.size();i++) {
	      model_diag_method.push_back(atoi(linedata[i].c_str()));
	    }
	    lineok=1;
	  }
	}
	if (linedata[0]==string("WriteHamiltToFile")) {
	  if (linedata.size()>=2) {
	    for (i=1;i<linedata.size();i++) {
	      model_hamilttofile.push_back(linedata[i]);
	    }
	    lineok=1;
	  }
	}
	if (linedata[0]==string("LineHamiltonian")) {
	  if (linedata.size()==4) {
	    model_linehamilt.push_back(atoi(linedata[1].c_str()));
	    model_linehamilt.push_back(atoi(linedata[2].c_str()));
	    model_linehamilt.push_back(atoi(linedata[3].c_str()));
	    lineok=1;
	  }
	}
	if (linedata[0]==string("PrintBasis")) {
	  if (linedata.size()>=2) { 
	    lineok=1;
	    for (i=1;i<linedata.size();i++) {
	      if ((linedata[i]==string("y")) || (linedata[i]==string("n"))) {
		model_printbasis.push_back(linedata[i]);
	      } else {
		printf("The flags for printing out the basis have to be ");
		printf("either 'y' or 'n'.\n");
		lineok=0;
	      }
	    }
	  }
	}
	if (linedata[0]==string("PrintEigenvalues")) {
	  if (linedata.size()>=2) { 
	    lineok=1;
	    for (i=1;i<linedata.size();i++) {
	      if ((linedata[i]==string("y")) || (linedata[i]==string("n")) ||
		  (linedata[i]==string("f"))) {
		model_printeigenvals.push_back(linedata[i]);
	      } else {
		printf("The flags for printing out the eigenvalues have to be ");
		printf("either 'y', 'n' or 'f'.\n");
		lineok=0;
	      }
	    }
	  }
	}
	if (linedata[0]==string("PrintEigenvectors")) {
	  if (linedata.size()>=2) { 
	    lineok=1;
	    for (i=1;i<linedata.size();i++) {
	      if ((linedata[i]==string("y")) || (linedata[i]==string("n")) ||
		  (linedata[i]==string("f"))) {
		model_printeigenvecs.push_back(linedata[i]);
	      } else {
		printf("The flags for printing out the eigenvectors have to be ");
		printf("either 'y', 'n' or 'f'.\n");
		lineok=0;
	      }
	    }
	  }
	}
	if (linedata[0]==string("PrintInputsummary")) {
	  if (linedata.size()==2) {
	    if ((linedata[1]==string("y")) || (linedata[1]==string("n"))) {
	      model_printinputsummary=linedata[1];
	      lineok=1;
	    }
	  }
	}
	if (linedata[0]==string("VerboseLevel")) {
	  if (linedata.size()==2) {
	    model_verboselevel=atoi(linedata[1].c_str());
	    lineok=1;
	  }
	}
      }

      if (!lineok==1) {
	printf("Error in inputfile in line %d! Aborting ..\n",linecount);
	exit(0);
      }
    }

    getline(fin,inputline,delimiter);
    linecount=linecount+1;
  }

  // do some post processing
  for (i=0;i<model_ptmaxnums.size();i++) {
    if (model_ptmaxnums[i]==-1) { model_ptmaxnums[i]=model_totalparticlenum; };
  }
  for (i=0;i<tmphnames.size();i++) {
    nameiter=hamilt_names.find(tmphnames[i]);
    if (nameiter==hamilt_names.end()) {
      printf("Unknown Hamiltonian : %s. Aborting ..\n",tmphnames[i].c_str());
      exit(0);
    } else {
      model_usehamiltonian.push_back(nameiter->second);
    }
  }

  check_parameter();

  fin.close();
}

void calc()
{
  int i,j,Un,N,J,P,maxpj,flag,isfmaxn;
  char fname[1000];
  CDidParticleCFP *cfptable;
  CLinearBasisQN2 bas1,bas2,bas3;
  CDReducedMatrix2 *rmats1,*rmats2,*rmats3,*rmats4;
  CDReducedMatrix2 rmat12,rmat34,rmat1234;
  vector<int> js1,js2,js3;
  vector<int> terms,tmp;
  vector<int> basis_limits2,basis_limits3;
  vector<double> coefs;
  vector<vector<int> > sensbas2,sensbas3;
  int pt1,pj1,pt2,pj2,pt3,pj3,pt4,pj4,interj,j1,j2,j3a,j3b,oj;
  int icalc,opindex,hindex;

  CDReducedMatrix2 hsparsemat,osparsemat;

  cfptable=new CDidParticleCFP[model_ptnum];
  // read the ISF's
  if (model_verboselevel>0) {
    printf("Reading the ISF's from harddisk ..\n");
  }
  for (i=0;i<model_ptnum;i++) {
    cfptable[i].Clear();
    Un=model_ptsymmetries[i];

    if (model_filterseniorities[i]==-1) {
      isfmaxn=model_ptmaxnums[i];
    } else {
      if (model_ptmaxnums[i]<model_filterseniorities[i]+2) {
	isfmaxn=model_ptmaxnums[i];
      } else {
	isfmaxn=model_filterseniorities[i]+2;
      }
    }
    //    isfmaxn=model_ptmaxnums[i];

    for (j=1;j<=isfmaxn;j++) { // j is the particle number
      if (!((Un==1) && (j>1))) { // s-boson states do only have nu=0 or nu=1
	//	sprintf(fname,"/zeta/heinze/ISFDataBase/disf_un%02d_n%02d.dat",Un,j);
	sprintf(fname,"%s/disf_un%02d_n%02d.dat",model_isf_dir.c_str(),Un,j);
	if (model_verboselevel>0) {
	  //	  printf("/zeta/heinze/ISFDataBase/disf_un%02d_n%02d.dat .. ",Un,j);
	  printf("Reading %s/disf_un%02d_n%02d.dat .. ",model_isf_dir.c_str(),Un,j);
	}
	cfptable[i].MergeDFile(fname);
	if (model_verboselevel>0) { printf("done.\n"); };
      }
    }
  }

  // get the highest spin of all distinguishable particles
  maxpj=0;
  for (i=0;i<model_ptsymmetries.size();i++) {
    if (model_ptsymmetries[i]-1>maxpj) { maxpj=model_ptsymmetries[i]-1; };
  }

  allhsparsemats.clear();
  allhsizes.clear();
  allbasis.clear();
  allobsparsemats.clear();
  allobsparsematsptr.clear();
  N=model_totalparticlenum;

  if (model_verboselevel>0) {
    printf("Starting to calculate all hamiltonian matrices.\n");
  }

  for (icalc=0;icalc<model_totalspin.size();icalc++) {
    if (model_verboselevel>0) {
      printf("Constructing hamiltonian matrix %d of %d :\n",
	     icalc+1,(int)model_totalspin.size());
    }
    J=model_totalspin[icalc];
    P=model_totalparity[icalc];

    hindex=model_usehamiltonian[icalc];

    if (model_verboselevel>0) { printf("J = %d/2 ; P =% d\n",J,P); };

    if (model_verboselevel>1) {
      printf("Starting to construct the basis 1 of 3 ..\n");
    }

    bas1.ConstructBasis(N,J,J,P,model_ptsymmetries,model_ptmaxnums,
			model_basislimits,model_ptparities,
			model_filterseniorities);

    // TODO : I can make the intermediate bases (bas2 and bas3)
    // much smaller. I can give a vector with particle type maxnums
    // which has smaller numbers. The maxnums can be reduced by 1 for bas2
    // and they can be reduced by 2 for bas3. This will bring me an
    // advantage of course only if the basis bas1 was restricted by some
    // maxnum which is smaller than N.

    if (J-maxpj<0) { j=0; } else { j=J-maxpj; };
    // Prepare the vectors which hold the allowed seniorities for the bas2
    sensbas2.clear();
    for (i=0;i<model_filterseniorities.size();i++) {
      tmp.clear();
      if (model_filterseniorities[i]!=-1) {
	tmp.push_back(model_filterseniorities[i]+1);
	tmp.push_back(model_filterseniorities[i]);
	if (model_filterseniorities[i]-1>=0) {
	  tmp.push_back(model_filterseniorities[i]-1);
	}
      }
      sensbas2.push_back(tmp);
    }

    if (model_verboselevel>1) {
      printf("Starting to construct the basis 2 of 3 ..\n");
    }

    // Eliminate the '=' basis restrictions. You want to understand why?
    // Think about it!
    basis_limits2.clear();
  /*
    for (i=0;i<model_basislimits.size();i++) {
      //      printf("%d ",model_basislimits[i]);
      if ((i+1)%(model_ptnum+2)==model_ptnum+1) {
	      if (model_basislimits[i]!=1) {
	         for (j=0;j<model_ptnum;j++) {
	            basis_limits2.push_back(model_basislimits[i-model_ptnum+j]);
	         }
	      basis_limits2.push_back(model_basislimits[i]);
	      basis_limits2.push_back(model_basislimits[i+1]);
	      }
      }
      //      if ((i+1)%(model_ptnum+2)==0) { printf("\n"); };
    }
    //    printf("%d\n",basis_limits2.size());
   */

    //    bas2.ConstructBasis(N-1,j,J+maxpj,0,model_ptsymmetries,model_ptmaxnums,
    //    			model_basislimits,model_ptparities);
    //    bas2.ConstructBasis(N-1,j,J+maxpj,0,model_ptsymmetries,model_ptmaxnums,
    //			model_basislimits,model_ptparities,sensbas2);
    bas2.ConstructBasis(N-1,j,J+maxpj,0,model_ptsymmetries,model_ptmaxnums,
			basis_limits2,model_ptparities,sensbas2);

    if (J-2*maxpj<0) { j=0; } else { j=J-2*maxpj; };
    // Prepare the vectors which hold the allowed seniorities for the bas3
    sensbas3.clear();
    for (i=0;i<model_filterseniorities.size();i++) {
      tmp.clear();
      if (model_filterseniorities[i]!=-1) {
	tmp.push_back(model_filterseniorities[i]+2);
	tmp.push_back(model_filterseniorities[i]+1);
	tmp.push_back(model_filterseniorities[i]);
	if (model_filterseniorities[i]-1>=0) {
	  tmp.push_back(model_filterseniorities[i]-1);
	  if (model_filterseniorities[i]-2>=0) {
	    tmp.push_back(model_filterseniorities[i]-2);
	  }
	}
      }
      sensbas3.push_back(tmp);
    }

    if (model_verboselevel>1) {
      printf("Starting to construct the basis 3 of 3 ..\n");
    }

    //    bas3.ConstructBasis(N-2,j,J+2*maxpj,0,model_ptsymmetries,model_ptmaxnums,
    //			model_basislimits,model_ptparities);
    //    bas3.ConstructBasis(N-2,j,J+2*maxpj,0,model_ptsymmetries,model_ptmaxnums,
    //    			model_basislimits,model_ptparities,sensbas3);
    bas3.ConstructBasis(N-2,j,J+2*maxpj,0,model_ptsymmetries,model_ptmaxnums,
    			basis_limits2,model_ptparities,sensbas3);

    if (model_verboselevel>0) {
      printf("Number of states : %d\n",bas1.GetStateNum());
    }
    //    printf("Number of states of bas2 : %d\n",bas2.GetStateNum());
    //    printf("Number of states of bas3 : %d\n",bas3.GetStateNum());

    if (model_printbasis[icalc]==string("y")) {
      printf("Basis of calculation %d :\n",icalc+1);
      bas1.Print();
    }
//      cout<<endl;bas2.Print();cout<<endl;bas3.Print();cout<<endl;
    allbasis.push_back(bas1);
    allhsizes.push_back(bas1.GetStateNum());

    rmats1=new CDReducedMatrix2[model_ptnum];
    rmats2=new CDReducedMatrix2[model_ptnum];
    rmats3=new CDReducedMatrix2[model_ptnum];
    rmats4=new CDReducedMatrix2[model_ptnum];

    bas1.GetJArray(&js1);
    bas2.GetJArray(&js2);
    bas3.GetJArray(&js3);
//      cout<<"ptnum: "<<model_ptnum<<endl;
    if (model_verboselevel>1) {
      printf("Contructing matrix representation of creation operators 1 of 2\n");
    }
    for (i=0;i<model_ptnum;i++) {
      rmats1[i].CalcReducedMatrix(&bas1,&bas2,&cfptable[i],i+1);
//        rmats1[i].Print();
    }
    if (model_verboselevel>1) {
      printf("Contructing matrix representation of creation operators 2 of 2\n");
    }
    for (i=0;i<model_ptnum;i++) {
      rmats2[i].CalcReducedMatrix(&bas2,&bas3,&cfptable[i],i+1);
    }
    if (model_verboselevel>1) {
      printf("Contructing matrix representation of annihilation operators 1 of 2\n");
    }
    for (i=0;i<model_ptnum;i++) {
      rmats3[i].CalcReducedMatrix(&bas3,&bas2,&cfptable[i],-(i+1));
    }
    if (model_verboselevel>1) {
      printf("Contructing matrix representation of annihilation operators 2 of 2\n");
    }
    for (i=0;i<model_ptnum;i++) {
      rmats4[i].CalcReducedMatrix(&bas2,&bas1,&cfptable[i],-(i+1));
//        rmats4[i].Print();
    }

    hsparsemat.SetNullMatrix(bas1.GetStateNum(),bas1.GetStateNum());

    hamiltonian[hindex].GetOneBodyTerms(&terms,&coefs);
    for (i=0;i<coefs.size();i++) {
      // IMPORTANT NOTE :
      // This works *ONLY* if we use for rmats1 and rmats4 the SAME pair of basis.
      pt1=terms[2*i+0];
      pj1=model_ptsymmetries[pt1-1]-1;
      pt2=terms[2*i+1];
      pj2=model_ptsymmetries[pt2-1]-1;
      rmat12.CalcReducedMatrix(&rmats1[pt1-1],&rmats4[pt2-1],&js1.front(),
			       &js2.front(),&js1.front(),pj1,pj2,0);
      rmat12.Mul3jAndPhase(J,-J,0,0,J,J,1);
      rmat12.Mul(coefs[i]);
      //      rmat12.Print();
      hsparsemat=hsparsemat+rmat12;
    }

    hamiltonian[hindex].GetTwoBodyTerms(&terms,&coefs);
    for (i=0;i<coefs.size();i++) {
      //    printf("%d %d %d %d %d %f\n",terms[5*i],terms[5*i+1],terms[5*i+2],
      //	   terms[5*i+3],terms[5*i+4],coefs[i]);
      pt1=terms[5*i+0];
      pj1=model_ptsymmetries[pt1-1]-1;
      pt2=terms[5*i+1];
      pj2=model_ptsymmetries[pt2-1]-1;
      pt3=terms[5*i+2];
      pj3=model_ptsymmetries[pt3-1]-1;
      pt4=terms[5*i+3];
      pj4=model_ptsymmetries[pt4-1]-1;
      interj=terms[5*i+4];

      rmat12.CalcReducedMatrix(&rmats1[pt1-1],&rmats2[pt2-1],&js1.front(),
			       &js2.front(),&js3.front(),pj1,pj2,interj);
      rmat34.CalcReducedMatrix(&rmats3[pt3-1],&rmats4[pt4-1],&js3.front(),
			       &js2.front(),&js1.front(),pj3,pj4,interj);
      rmat1234.CalcReducedMatrix(&rmat12,&rmat34,&js1.front(),&js3.front(),
				 &js1.front(),interj,interj,0);

      rmat1234.Mul3jAndPhase(J,-J,0,0,J,J,1);
      rmat1234.Mul(coefs[i]);
      //      rmat1234.Print();
      hsparsemat=hsparsemat+rmat1234;
    }

    //    hsparsemat.Print();
    //    exit(0);

    if (model_hamilttofile.size()>0) {
      char *tmpname;
      printf("Writing hamiltonian to file : %s\n",model_hamilttofile[icalc].c_str());
      tmpname=new char[model_hamilttofile[icalc].size()+1];
      sprintf(tmpname,"%s",model_hamilttofile[icalc].c_str());
      hsparsemat.WriteToFile(tmpname);
      delete [] tmpname;
    }

    allhsparsemats.push_back(hsparsemat);

    // clean up
    delete [] rmats1;
    delete [] rmats2;
    delete [] rmats3;
    delete [] rmats4;
  }

  if (model_hamilttofile.size()>0) {
    printf("The programm was told to write the hamiltonian matrices to files,\n");
    printf("which was done. It is assumed that these files are needed for\n");
    printf("something this program is not able to do, thus it exits here and\n");
    printf("wishes the user good luck.\n");
    exit(0);
  }

  if (model_verboselevel>0) {
    printf("Finished to calculate all hamiltonian matrices.\n");
  }

  if (model_verboselevel>0) {
    printf("Starting to calculate matrices of operators ..\n");
  }
  // Now we calculate the matrices for the other operators (one body)
  if (model_linehamilt.size()==3) { // if we have to do a line hamiltonian
    if (matels_to_calc.size()/5>0) { // and have some other operatos
      printf("The calculation of operator exspectation values together\n");
      printf("with a line hamiltonian is not implemented yet! Aborting ..\n");
      exit(0);
    }
  }

  for (icalc=0;icalc<matels_to_calc.size()/5;icalc++) {
    if (model_verboselevel>0) {
      printf("Constructing operator matrix %d of %d ..\n",
	     icalc+1,(int)matels_to_calc.size()/5);
    }

    // check if we already calculated the appropriate operator
    // matrix and if yes, set only the pointer to it
    flag=-1;
    for (i=0;i<icalc;i++) {
      if (flag==-1) {
	if ((matels_to_calc[icalc*5]==matels_to_calc[i*5]) &&
	    (matels_to_calc[icalc*5+2]==matels_to_calc[i*5+2]) &&
	    (matels_to_calc[icalc*5+3]==matels_to_calc[i*5+3])) {
	  //	  flag=i;
	  flag=allobsparsematsptr[i];
	}
      }
    }
    // If we calculated the appropriate operator matrix already, 'flag' is
    // the correct index, flag==-1 otherwise.

    if (flag==-1) {
      bas1=allbasis[matels_to_calc[icalc*5]-1];
      bas2=allbasis[matels_to_calc[icalc*5+3]-1];
      osparsemat.SetNullMatrix(bas1.GetStateNum(),bas2.GetStateNum());
      // Calculate the appropriate intermediate basis
      // Get the spin of the basis
      j1=model_totalspin[matels_to_calc[icalc*5]-1];
      j2=model_totalspin[matels_to_calc[icalc*5+3]-1];
      // Calculate the highest and lowest spins for the intermediate basis
      if ((j1-maxpj)>=(j2-maxpj)) {
	j3a=j1-maxpj;
      } else {
	j3a=j2-maxpj;
      }
      if (j3a<0) { j3a=0; };
      if ((j1+maxpj)<=(j2+maxpj)) {
	j3b=j1+maxpj;
      } else {
	j3b=j2+maxpj;
      }

      //      printf("test1\n");
      /*
      printf("%d %d %d\n",N-1,j3a,j3b);
      for (i=0;i<model_ptsymmetries.size();i++) {
	printf("%d ",model_ptsymmetries[i]);
      }
      printf("\n");
      for (i=0;i<model_ptmaxnums.size();i++) {
	printf("%d ",model_ptmaxnums[i]);
      }
      printf("\n");
      for (i=0;i<model_basislimits.size();i++) {
	printf("%d ",model_basislimits[i]);
      }
      printf("\n");
      for (i=0;i<model_ptparities.size();i++) {
	printf("%d ",model_ptparities[i]);
      }
      printf("\n");
      */
      //      model_basislimits[4]=2;
      //      model_basislimits[9]=2;

      basis_limits3.clear();

      // TODO : filter seniorities !!!!!!
	//      bas3o.ConstructBasis(N-1,j3a,j3b,0,model_ptsymmetries,model_ptmaxnums,
	//			  model_basislimits,model_ptparities);
      bas3.ConstructBasis(N-1,j3a,j3b,0,model_ptsymmetries,model_ptmaxnums,
			  basis_limits3,model_ptparities);

      //      printf("test2\n");
      /*
      printf("test\n");
      bas1.Print();
      printf("test\n");
      bas2.Print();
      printf("test\n");
      bas3.Print();
      exit(0);
      */

      bas1.GetJArray(&js1);
      bas2.GetJArray(&js2);
      bas3.GetJArray(&js3);

      // Calculate the reduced matrix element of creators and annihilators
      rmats1=new CDReducedMatrix2[model_ptnum];
      rmats2=new CDReducedMatrix2[model_ptnum];

      for (i=0;i<model_ptnum;i++) {
	rmats1[i].CalcReducedMatrix(&bas1,&bas3,&cfptable[i],i+1);
	rmats2[i].CalcReducedMatrix(&bas3,&bas2,&cfptable[i],-(i+1));
      }

      opindex=matels_to_calc[icalc*5+2];

      for (i=0;i<allobopscoefs[opindex].size();i++) { // do term by term
	pt1=allobops[opindex][i*3];
	pt2=allobops[opindex][i*3+1];
	pj1=model_ptsymmetries[pt1-1]-1;
	pj2=model_ptsymmetries[pt2-1]-1;
	oj=allobops[opindex][i*3+2];
	//      printf("coef = %e\n",allobopscoefs[opindex][i]);
	//      printf("%d %d %d %d %d\n",pt1,pt2,pj1,pj2,oj);
	rmat12.CalcReducedMatrix(&rmats1[pt1-1],&rmats2[pt2-1],&js1.front(),
				 &js3.front(),&js2.front(),pj1,pj2,oj);
	//	rmat12.Print();
	rmat12.Mul(allobopscoefs[opindex][i]);
	//	rmat12.Mul3jAndPhase(j1,-j1,oj,j1-j2,j2,j2,1); ????
	osparsemat=osparsemat+rmat12;
      }

      allobsparsemats.push_back(osparsemat);
      //      osparsemat.Print();
      //      allobsparsematsptr.push_back(icalc); error!
      { // for testing write the operator matrix to a file
	//	printf("Writing operator matrix to file : test.txt\n");
	//	osparsemat.WriteToFile("test.txt");
      }

      allobsparsematsptr.push_back(allobsparsemats.size()-1);

      delete [] rmats1;
      delete [] rmats2;

    } else { // we already calculated the matrix, just set the index
      allobsparsematsptr.push_back(flag);
    }
  }

  if (model_verboselevel>0) {
    printf("Finished to calculate all operator matrices.\n");
  }

  //  for (i=0;i<allobsparsematsptr.size();i++) {
  //    printf("%d , ",allobsparsematsptr[i]);
  //  }
  //  printf("\n");

  delete [] cfptable;
}

void calc_eigenvals()
{
  int i,j,k,hsize,icalc;
  char buf[1024];
  CSubroutines subs;
  vector<double> hdensemat,eigenvals,eigenvecs;

  // for the line hamiltonian
  int pnum,statenum;
  double lambda,energy0;
  CDReducedMatrix2 lhamilt1,lhamilt2,lhamilt;
  vector<double> lineres;

  alleigenvals.clear();
  alleigenvecs.clear();

  if (model_linehamilt.size()==0) { // no line hamiltonian, just normal calculations
    if (model_verboselevel>0) {
      printf("Doing the normal calculations ..\n");
    }
    for (icalc=0;icalc<allhsparsemats.size();icalc++) {

      if (allhsparsemats[icalc].GetRowNum()!=allhsparsemats[icalc].GetColumnNum()) {
	printf("Hamiltonian is not a square matrix!?!? Aborting ..\n");
	exit(0);
      }

      hsize=allhsparsemats[icalc].GetRowNum();

      // just to test something
      // printf("Size = %d x %d = %d ; NonZeros = %d ; Sparsness = %f\n",hsize,hsize,
      //	   hsize*hsize,allhsparsemats[icalc].NonZeros,
      //	   1.0-(double)(allhsparsemats[icalc].NonZeros)/(double)(hsize*hsize));

      if (model_states_to_calc[icalc]>hsize) {
	model_states_to_calc[icalc]=hsize;
      }

      if (model_diag_method[icalc]==1) {
	// full diagonalisation with the built in Jacobi procedure
	hdensemat.resize(hsize*hsize);
	//    sprintf(buf,"hsparsematJ%dP%d.txt",
	//	    model_totalspin[icalc],model_totalparity[icalc]);
	//    allhsparsemats[icalc].WriteToFile(buf);

	allhsparsemats[icalc].GetMatrix(&hdensemat.front());

	// print out the hamilton matrix
	//	for (i=0;i<hsize;i++) {
	//	  for (j=0;j<hsize;j++) { printf("%f ",hdensemat[i*hsize+j]); };
	//	  printf("\n");
	//	}

	subs.SymmetricDiag(&hdensemat,&eigenvals,&eigenvecs);
    
	if (model_states_to_calc[icalc]!=-1) { // store not all eigenvals and eigenvecs
	  eigenvals.erase(eigenvals.begin()+model_states_to_calc[icalc],
			  eigenvals.end());
	  eigenvecs.erase(eigenvecs.begin()+hsize*model_states_to_calc[icalc],
			  eigenvecs.end());
	}
	alleigenvals.push_back(eigenvals);
	alleigenvecs.push_back(eigenvecs);
      }

   /*
   if (model_diag_method[icalc]==2) {
	   CArpackWrapper arpack;
	   int eigenvalnum;
	   eigenvalnum=model_states_to_calc[icalc];
	   if (eigenvalnum==-1) { eigenvalnum=hsize; };
	   // prepare space for the eigenvalues and eigenvectors
	   eigenvals.resize(eigenvalnum);
	   eigenvecs.resize(eigenvalnum*hsize);
	   // if we need all eigenvalues, split the calculation into two parts
	   if (eigenvalnum==hsize) { // we want ALL eigenvalus
	     arpack.Arpack_dsaupd(hsize/2,&eigenvals.front(),&eigenvecs.front(),1,0,
			       allhsparsemats[icalc].NCols,
			       allhsparsemats[icalc].IMatrix,
			       allhsparsemats[icalc].JMatrix,
			       allhsparsemats[icalc].Matrix);
	   arpack.Arpack_dsaupd(hsize-hsize/2,&eigenvals.front()+hsize/2,
			       &eigenvecs.front()+hsize*(hsize/2),1,1,
			       allhsparsemats[icalc].NCols,
			       allhsparsemats[icalc].IMatrix,
			       allhsparsemats[icalc].JMatrix,
			       allhsparsemats[icalc].Matrix);
	   } else { // not all eigenvalues
	      arpack.Arpack_dsaupd(eigenvalnum,
			       &eigenvals.front(),&eigenvecs.front(),1,0,
			       allhsparsemats[icalc].NCols,
			       allhsparsemats[icalc].IMatrix,
			       allhsparsemats[icalc].JMatrix,
			       allhsparsemats[icalc].Matrix);
	   }
	alleigenvals.push_back(eigenvals);
	alleigenvecs.push_back(eigenvecs);

	//	printf("This is not implemented yet! Aborting ..\n");
	//	exit(0);
   }
   */

      if (model_printeigenvals[icalc]==string("y")) {
	printf("Eigenvalues for J = %d/2 ; P = %d ; (%d of %d) :\n",
	       model_totalspin[icalc],model_totalparity[icalc],
	       (int)alleigenvals[icalc].size(),hsize);
	for (i=0;i<alleigenvals[icalc].size();i++) {
	  // substract the groundstate
	  //	  printf("%f ",alleigenvals[icalc][i]-alleigenvals[0][0]);
	  if (model_groundstateindex==-1) {
	    energy0=0;
	  } else {
	    energy0=alleigenvals[0][0];
	  }
	  printf("%f ",alleigenvals[icalc][i]-energy0);
	  //	  if ((i+1)%6==0) { printf("\n"); };
	}
	printf("\n");
      }
      if (model_printeigenvals[icalc]==string("f")) {
	for (i=0;i<alleigenvals[icalc].size();i++) {
	  printf("EVAL %d %d %d %d %d %d %.20f\n",icalc+1,model_totalspin[icalc],
		 model_totalparity[icalc],i+1,(int)alleigenvals[icalc].size(),
		 hsize,alleigenvals[icalc][i]);
	}
	printf("\n");
      }

      if (model_printeigenvecs[icalc]==string("y")) {
	printf("Eigenvectors for J = %d/2 ; P = %d ; (%d of %d) :\n",
	       model_totalspin[icalc],model_totalparity[icalc],
	       (int)alleigenvals[icalc].size(),hsize);
	for (i=0;i<alleigenvals[icalc].size();i++) {
	  printf("|%d> = ",i);
	  for (j=0;j<hsize;j++) {
	    printf("%e ",alleigenvecs[icalc][i*hsize+j]);
	  }
	  printf("\n");
	}
	printf("\n");
      }
      if (model_printeigenvecs[icalc]==string("f")) {
	for (i=0;i<alleigenvals[icalc].size();i++) {
	  for (j=0;j<hsize;j++) {
	    printf("EVEC %d %d %d %d %d %d %d %e\n",icalc+1,model_totalspin[icalc],
		   model_totalparity[icalc],i+1,j+1,(int)alleigenvals[icalc].size(),
		   hsize,alleigenvecs[icalc][i*hsize+j]);
	  }
	}
	printf("\n");
      }

    }
  }

  if (model_linehamilt.size()==3) {
    printf("#Doing a line hamiltonian :\n");
    lineres.clear();
    lhamilt1=allhsparsemats[model_linehamilt[0]-1];
    lhamilt2=allhsparsemats[model_linehamilt[1]-1];
    pnum=model_linehamilt[2];
    hsize=lhamilt1.GetRowNum();
    statenum=model_states_to_calc[model_linehamilt[0]-1];
    if (statenum==-1) { statenum=hsize; };
    if (statenum>hsize) { statenum=hsize; };
    for (icalc=0;icalc<pnum;icalc++) {
      lambda=(double)(icalc)/(double)(pnum-1);
      lhamilt=lhamilt1*(1.0-lambda)+lhamilt2*lambda;
      printf("#%f\n",lambda);
      // If the standart output is piped to somewhere, we still see
      // the status on the shell with the following line if standart error
      // is not also redirected.
      //      cerr << lambda << endl;

      if (lhamilt.GetRowNum()!=lhamilt.GetColumnNum()) {
	printf("Hamiltonian is not a square matrix!?!? Aborting ..\n");
	exit(0);
      }

      if (model_diag_method[0]==1) {
	hdensemat.resize(hsize*hsize);
	lhamilt.GetMatrix(&hdensemat.front());
	subs.SymmetricDiag(&hdensemat,&eigenvals,&eigenvecs);
      }
   /*
      if (model_diag_method[0]==2) {
	CArpackWrapper arpack;
	eigenvals.resize(hsize);
	eigenvecs.resize(1); // just to have a valid pointer,we dont get eigenvectors
	arpack.Arpack_dsaupd(hsize/2,&eigenvals.front(),&eigenvecs.front(),0,0,
			     lhamilt.NCols,lhamilt.IMatrix,lhamilt.JMatrix,
			     lhamilt.Matrix);
	arpack.Arpack_dsaupd(hsize-hsize/2,&eigenvals.front()+hsize/2,
			     &eigenvecs.front(),0,1,
			     lhamilt.NCols,lhamilt.IMatrix,lhamilt.JMatrix,
			     lhamilt.Matrix);	
      }
   */

      // printf("%e\n",eigenvals[0]);
      lineres.push_back(lambda);
      for (i=0;i<statenum;i++) { lineres.push_back(eigenvals[i]); };
    }

    // write output
    for (j=1;j<=statenum;j++) {
      for (i=0;i<lineres.size()/(statenum+1);i++) {
	printf("%e %e\n",lineres[i*(statenum+1)+0],lineres[i*(statenum+1)+j]);
	// lineres[i*(statenum+1)+j]-lineres[i*(statenum+1)+1]);
      }
      printf("\n");
    }

    // do the finite temperature thing (idea of Pavel Cejnar)
    /*
    {
      int it;
      for (it=1;it<=1;it++) {
	int pointnum=lineres.size()/(statenum+1);
	vector<double> partfuncs(pointnum,0.0); // to store the partition functions
	vector<double> ufuncs(pointnum*statenum,0.0); // to store U_i(\lambda)
	vector<double> cfuncs(pointnum*statenum,0.0); // to store C_i(\lambda)
	vector<double> finalcs(pointnum,0.0);
	double z,u,c,c1,c2,t,lambdadiff;
	t=it*1.0; // choose a proper value for the temperature
	// calculate all partition functions
	for (j=0;j<pointnum;j++) { // j is the point index
	  z=0;
	  for (i=1;i<=statenum;i++) { z=z+exp(-lineres[j*(statenum+1)+i]/t); };
	  partfuncs[j]=z;
	  //	printf("partiction function : %d %e\n",j,z);
	}
	// calculate all U_i(\lambda)
	for (j=0;j<pointnum;j++) { // j is the point index
	  for (i=1;i<=statenum;i++) {
	    // calculate U_i(\lambda) where lambda is given by the index j
	    u=0;
	    for (k=1;k<=statenum;k++) {
	      if (k!=i) {
		u=u-log(fabs(lineres[j*(statenum+1)+i]-lineres[j*(statenum+1)+k]));
		// printf("%d %d %d %e %e %e %e\n",j,j,k,lineres[j*(statenum+1)+i],
		//       lineres[j*(statenum+1)+k],
		//       fabs(lineres[j*(statenum+1)+i]-lineres[j*(statenum+1)+k]),u);
	      }
	    }
	    ufuncs[j*statenum+i-1]=u;
	    //	    printf("u : %d %f %d %e\n",j,lineres[j*(statenum+1)],i,u);
	  }
	}
	// calculate all C_i(\lambda)
	lambdadiff=fabs(lineres[0]-lineres[statenum+1]);
	for (j=2;j<pointnum-2;j++) { // j is the point index
	  for (i=1;i<=statenum;i++) {
	    // calculate C_i(\lambda)

	    //c1=(lineres[j*(statenum+1)]*ufuncs[j*statenum+i-1]-
	    //	lineres[(j-2)*(statenum+1)]*ufuncs[(j-2)*statenum+i-1])/lambdadiff;
	    //c2=(lineres[(j+2)*(statenum+1)]*ufuncs[(j+2)*statenum+i-1]-
	    //	lineres[j*(statenum+1)]*ufuncs[j*statenum+i-1])/lambdadiff;
	    c1=(ufuncs[j*statenum+i-1]-ufuncs[(j-2)*statenum+i-1])/lambdadiff;
	    c2=(ufuncs[(j+2)*statenum+i-1]-ufuncs[j*statenum+i-1])/lambdadiff;
	    //	    c=-1.0/(statenum-1)*lineres[j*(statenum+1)]*(c2-c1)/lambdadiff;
	    c=-1.0/(statenum-1)*(c2-c1)/lambdadiff;
	    cfuncs[j*statenum+i-1]=c;
	    // printf("c : %d %d %e %e %e %e %e\n",j,i,c1,c2,c,ufuncs[j*statenum+i-1],
	    //	      ufuncs[(j-2)*statenum+i-1]);
	  }
	}
	// calculate C(\lambda)
	for (j=2;j<pointnum-2;j++) {
	  c=0;
	  for (i=0;i<statenum;i++) {
	    c=c+exp(-lineres[j*(statenum+1)+i+1]/t)*cfuncs[j*statenum+i+1];
	  }
	  finalcs[j]=c/partfuncs[j];
	}
	// output

	//	for (j=2;j<pointnum-2;j++) {
	//	  printf("%e %e %e\n",lineres[j*(statenum+1)],t,finalcs[j]);
	//	}
	//	printf("\n");

	//	for (i=0;i<statenum;i++) {
	//	  for (j=2;j<pointnum-2;j++) {
	//	    printf("%e %e %d %e\n",t,lineres[j*(statenum+1)],i,
	//		   exp(-lineres[j*(statenum+1)+1+i]/t)/partfuncs[j]);
	//	  }
	//	  printf("\n");

	for (i=0;i<statenum;i++) {
	  for (j=2;j<pointnum-2;j++) {
	    printf("%e %e %d %e\n",t,lineres[j*(statenum+1)],i,cfuncs[j*statenum+i]);
	  }
	  printf("\n");

	}
      }
    } // end finite temperature thing
    */

  }
}

void calc_op_matels()
{
  int imatel,i,size1,size2;
  double *outvec,*invec,*bravec;
  double matel;

  allobmatels.clear();

  for (imatel=0;imatel<matels_to_calc.size()/5;imatel++) {
    //    allobsparsemats[imatel].Print();

    size1=allhsizes[matels_to_calc[imatel*5]-1];
    size2=allhsizes[matels_to_calc[imatel*5+3]-1];

    invec=&alleigenvecs[matels_to_calc[imatel*5+3]-1].front()+
      size2*(matels_to_calc[imatel*5+4]-1);

    outvec=new double[size1];

    //    for (i=0;i<size2;i++) {
    //      printf("%f ",invec[i]);
    //    }
    //    printf("\n");

    allobsparsemats[allobsparsematsptr[imatel]].MulVec(invec,outvec);

    //    for (i=0;i<size1;i++) {
    //      printf("%e ",outvec[i]);
    //    }
    //    printf("\n");

    // Scalar product with the bra-vector
    //    bravec=&alleigenvecs[matels_to_calc[imatel*5+0]-1].front()+
    //      size2*(matels_to_calc[imatel*5+1]-1);
    bravec=&alleigenvecs[matels_to_calc[imatel*5+0]-1].front()+
      size1*(matels_to_calc[imatel*5+1]-1);
    matel=0;
    //    for (i=0;i<size1;i++) {
    //      printf("%e ",bravec[i]);
    //    }
    //    printf("\n");
    for (i=0;i<size1;i++) {
      matel=matel+outvec[i]*bravec[i];
    }

    printf("<%d %d||%s||%d %d>=%e  ; ||^2=%e ; ||^2/(2J_i+1)=%e\n",
	   matels_to_calc[imatel*5],
	   matels_to_calc[imatel*5+1],
	   allobopsnames_array[matels_to_calc[imatel*5+2]].c_str(),
	   matels_to_calc[imatel*5+3],matels_to_calc[imatel*5+4],
	   matel,matel*matel,
	   matel*matel/(model_totalspin[matels_to_calc[imatel*5]-1]+1));

    allobmatels.push_back(matel);

    delete [] outvec;
  }
}

int main(int argc,char *argv[])
{
  if (argc!=2) {
    printf("Wrong number of arguments!\n");
    printf("You have to give an input file. Aborting ..\n");
    exit(0);
  }

  read_inputfile(argv[1]);

  if (model_printinputsummary==string("y")) {
    print_input_summary();
  }

  calc();

  calc_eigenvals();

  calc_op_matels();
}
