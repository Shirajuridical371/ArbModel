# ArbModel

Modified version of the **ArbModel** nuclear structure code originally developed by  
**S. Heinze**, *University of Cologne (2008)*.

> Citation:  
> S. Heinze, *computer program ArbModel*, University of Cologne (2008).

All original source files are included with minimal changes, except for a major optimization of the
`SymmetricDiag` routine in `csubroutines.cpp`, which now uses the [Eigen](https://eigen.tuxfamily.org)
C++ linear algebra library for symmetric matrix diagonalization.

---

## ✨ Key Improvement

The modified `SymmetricDiag` function employs Eigen’s `SelfAdjointEigenSolver`, resulting in a
**~100× performance improvement** compared to the legacy LAPACK-based method on modern CPUs.

- **Tested on:** Apple M3 MacBook Pro  
- **Compiler:** `clang++` (Apple Clang 16) with `-O3 -std=c++17`  
- **Dependency:** install Eigen via Homebrew  
  ```bash
  brew install eigen
  ```


### 🧠 Technical Note

The Eigen-based version preserves the same API as the original routine:

```cpp
void CSubroutines::SymmetricDiag(
    std::vector<double>* mat,
    std::vector<double>* eigenvals,
    std::vector<double>* eigenvecs)
{
    using namespace std;
    using namespace std::chrono;

    const bool verbose = (std::getenv("ARB_VERBOSE") && string(std::getenv("ARB_VERBOSE")) == "1");

    const size_t nn = mat->size();
    const int N = (int)llround(std::sqrt((double)nn));
    if ((size_t)N * (size_t)N != nn) {
        cerr << "CSubroutines::SymmetricDiag(..): matrix not square!\n";
        std::exit(1);
    }

    time_point<high_resolution_clock> t0;
    if (verbose) {
        cout << "SymmetricDiag: matrix dimension = " << N
             << " (" << (size_t)N * (size_t)N << " elements, ≈ "
             << ((double)N * (double)N * 8.0 / 1e6) << " MB)\n";
        t0 = high_resolution_clock::now();
    }

    // Map the row-major flat buffer into an Eigen matrix (no copy)
    Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>
        A(mat->data(), N, N);

    // Solve: self-adjoint (symmetric) eigendecomposition
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(A);
    if (es.info() != Eigen::Success) {
        cerr << "Eigen decomposition failed\n";
        std::exit(2);
    }

    // Copy back results
    const auto& vals = es.eigenvalues();    // size N, ascending
    const auto& vecs = es.eigenvectors();   // NxN, columns = eigenvectors

    eigenvals->assign(vals.data(), vals.data() + N);
    eigenvecs->resize((size_t)N * (size_t)N);

    // Flatten eigenvectors (column j) into your row-major 1D layout
    for (int j = 0; j < N; ++j)
        for (int i = 0; i < N; ++i)
            (*eigenvecs)[(size_t)j * (size_t)N + (size_t)i] = vecs(i, j);

    if (verbose) {
        auto t1 = high_resolution_clock::now();
        auto ms = duration_cast<milliseconds>(t1 - t0).count();
        cout << "SymmetricDiag (Eigen) completed in "
             << ms << " ms (" << (ms / 1000.0) << " s)\n";
    }
}
```

### 🧩 Original Implementation (Commented Out)

The original ArbModel diagonalization code remains in the source, commented out for reference and validation.

### 🧮 Performance Note

Benchmarks show Eigen’s `SelfAdjointEigenSolver` achieves over **100× speedup** for typical large
symmetric matrices (N ≥ 1000) compared to the original solver used in ARBMODEL.

---

## ⚙️ Build Instructions

Use the provided `compile.sh`. Modification may be needed for your machine.

---

## ⚖️ License and Attribution

This repository is distributed for **academic and research use** under the same spirit as the original **ARBMODEL** code.

Please cite:

> S. Heinze, *computer program ARBMODEL*, University of Cologne (2008).

Additional acknowledgment is **welcome but not required**:

> Modifications by A. L. Conley (https://github.com/alconley), implementing Eigen-based diagonalization for high-performance computations (2025).

---

### 🌐 About This Project

This repository to make **nuclear physics tools** more accessible and to reproduce/obtain easier.

---

# Manual for the Numerical Code (ArbModel)

This part is copied from Dr. Heinze dissertation.

---

## A.1 What the code can do

The code ArbModel can calculate eigenvalues and eigenvectors as well as matrix elements of
operators according to the previously calculated eigenstates of a systems of identical and dis-
tinguishable particles (bosons and/or fermions) with a given spin. The number of particles can
be choosen freely. The most general hamiltonian given in the formalism of second quantization
including one and two body terms can be diagonalized. The user is bound to seniorities of iden-
tical particles for which a table of isoscalar factors exists on the harddisk. Another restriction
is simply the calculation time which depends strongly on the size of the hamiltonian matrix.

---

## A.2 A simple example: sd-IBM1 calculations

In this section a simple example for a sd-IBM1 calculation is given. This model is choosen
because of its simplicity. It consists only of two distinguishable particles with positive parity
and spin zero and two. The general hamiltonian has only six independent terms. An input file
for a calculation of 0+- and 2+-states including the B(E2; 2+_1 → 0+_1 ) will now be discussed in
detail line by line. This input file is a simple text file and should be given as an argument to
the program when it’s started from a command line. The first lines of the input file are:

```text
# tell the code where to find the isoscalar factors
FolderOfISF /somefolder/anotherfolder
```

The first line starts with a #, which does nothing else than starting a comment line. The user
is free to put as many comments as he wants at any position of the file as long as the comment
starts with the first character of the line. The parameter FolderOfISF is a string which is
the pathname of the folder where the isoscalar factors can be found. There is a special name
convention for the file containing isoscalar factors. For example “disf_un05_n07.dat” stands for:
double precision isoscalar factors for partciles with symmetry U(5) (d-bosons) and states with
seniority seven. In the present example the highest seniority for s-bosons is one and for d-bosons
ten. Thus there have to be the files disf_un05_n01.dat .. disf_un05_n10.dat and disf un01 n01.dat.
The specification of the model starts with the number of distinguishable particles.

```text
# the model:
ParticleTypeNumber 2
```

For spdf-IBM1 or sd-IBM2 this number would be four and for sdg-IBM1 there would be three
tpyes of particles. The two particle types of the sd-IBM1 are usually called s and d.

```text
ParticleTypeNames s d
```

These names are completely arbitrary strings of arbitrary length. The particle names will be
used later to give the hamiltonian and other operators explicitly in terms of creation and
annihilatoin operators. The program needs to know the parities as well as the spins of the
particles.

```text
ParticleTypeParities 1 1
ParticleTypeSymmetries 1 5
```

The parity can be either 1 or −1. The spin is given by the number of projections, i.e. the number
of magnetic substates. An odd number corresponds always to a boson and an even number to
a fermion always.
In the next two lines the restrictions of the basis are defined.

```text
ParticleTypeMaxNums -1 -1
ParticleTypeSeniorities -1 -1
```

A −1 means that there is no restriction. For an sdg-IBM1 calculation it can be useful to restrict
the number of g-bosons. If the hamiltonian conservs seniority for one or more particle types
a restriction of the seniority will result in a much smaller dimension of the hamiltonian which
will decrease the calculation time significantly.
In this example the user is interested in 0+ and 2+ states with total particle number of 10.
Thus the code has to calculate two hamiltonian matrices (one for spin zero and one for spin
two). Both sets of states have positive parity.

```text
# the states to calculate:
TotalParticleNumber 10
TotalSpin 0 4
TotalParity 1 1
```

There is only one total particle number for all calculations. The angular momenta are given
in units of hbar/2. The number of spins determines the number of calculations and the number
of parities has to be the same. If the dimension of the hamiltonian is large, one is usually
intersted in the first low lying states of each spin. Because of this there is the parameter
NumberOfStatesToCalc. A -1 means, that all eigenvalues will be calculated.

```text
NumberOfStatesToCalc -1 -1
```

It is possible to use a different hamiltonian for each calculation.

```text
UseHamiltonian H1 H1
```

The names of the hamiltonians are strings which can be choosen freely. The explicit form of these
hamiltonians has to be given later. Several methods are available to perform the diagonalization.
The parameter DiagMethod can be used to specify one method for each calculation.

```text
# tell the code which diagonalisation it should use
DiagMethod 1 1
```

There are several parameters which control the output on the screen. The parameters are more
or less self-explanatory except VerboseLevel. If it is set to a value bigger than zero, the user will
get additional information during the calculation process which is helpful only for real experts.

```text
# tell the code what to print out
PrintBasis y y
PrintEigenvalues y y
PrintEigenvectors n n
PrintInputsummary n
VerboseLevel 0
```

In this example the user is interested in the B(E2; 2+_1 → 0+_1 ) value. A definition of an operator
is done by putting an O as the first character of a line.

```text
# definition of some operatos
O E2 s + d ~ 4 1.0 d + s ~ 4 1.0
```

The second string is just a freely choosen name for the operator. The last line defines the E2-
operator as T(E2) := 1.0 ·(s† ×  ̃d)(4/2) + 1.0 ·(d† × ̃s)(4/2). Again all angular momenta are given
in units of hbar/2. There is no limit in the number of terms. The next task is to specify which
matrix elements of the previously defined operator which will be calculated.

```text
M 2 1 E2 1 1
```

The character M at the beginning of a line stands for “matrix element to calculate”. The first
two numbers determine the bra-state, the last two numbers the ket-state. A state is specified
by giving the number of the calculation and the number of the eigenstate of that calculation.
In this example the first calculation is for 0+ states and the second one for 2+ states. Thus
2 1 means the 2+_1 state and 1 1 means the 0+_1 state. Between the two states the name of the
operator is given. The output will be the reduced matrix element which is defined like in [10].
The last thing to do is to define the hamiltonian. This is done by giving one term in each line
(one-body or two-body). These terms are added internally to get the final hamiltonian operator.
The first character of such a line has to be a H followed by the name of the hamiltonian. Theses
names refer to the ones which were given by the parameter UseHamiltonian before.

```text
# the hamiltonian
H H1 d + d ~ 0.5
H H1 s + d ~ s + d ~ 4 -2.236e-01
H H1 s + d ~ d + s ~ 4 -2.236e-01
H H1 d + s ~ s + d ~ 4 -2.236e-01
H H1 d + s ~ d + s ~ 4 -2.236e-01
```

This defines the hamilton operator as

H := 0.5(d† ×  ̃d)(0) −0.2236[(s† ×  ̃d)^(4/2) ×(s† ×  ̃d)^(4/2)]^(0)−
0.2236[(s† ×  ̃d)^(4/2) ×(d† × ̃s)^(4/2)]^(0) − 0.2236[(d† × ̃s)^(4/2) ×(s† ×  ̃d)^(4/2)]^(0) − 0.2236[(d† × ̃s)^(4/2) ×(d† × ̃s)^(4/2)]^(0)]

The user has to ensure by his own that all operators are hermitian. The output of this example
would be the two bases for the 0+ and for the 2+ calculation and the two sets of eigenvalues
(energies). The code gives the absolute energies always, i.e. the pure eigenvalues from the
diagonalization.

---

## A.3 Input file reference

Every line of the input file has to be either empty or to start with one of the strings which will
be explained in detail now. There must be always least one blank character (space) after one of
these command strings. The user may add as many space characters (or tabulators) as he likes
after the parameter name or between the arguments.

### `#`
The whole line is treated as a **comment**.

### `FolderOfISF`
Gives the folder where the tables of isoscalar factors can be found.
**Example:**
```text
FolderOfISF /somefolder/anotherfolder
```

### `ParticleTypeNumber`
The argument to this parameter is one integer bigger than zero which gives the number
of distinguishable particles of the model.
**Example:**
```text
ParticleTypeNumber 2
```

### `ParticleTypeNames`
The parameter ParticleTypeNumber must be given before the code excepts this parameter.
The user has to give as many freely choosen strings here as there are distinguishable
particles of the model. These strings will be treated as names for the different particle
types. They can be used to give the hamiltonian as well as other operators. The strings
for the names can be of any length.
**Example:**
```text
ParticleTypeNames a b particle3 proton
```

### `ParticleTypeParities`
As ParticleTypeNames but here integers of either 1 or −1 have to be given.
**Example:**
```text
ParticleTypeParities 1 -1 1 -1
```

### `ParticleTypeSymmetries`
The arguments of this parameter determines the spins of all particles of the model. As
the name of the parameter suggests, the spin is given indirectly by the symmetry of the
particle. A particle with spin j (integer or half integer) is related to the symmetry U(n)
with n = 2j + 1. The arguments are exactly the dimension of the symmetry groups which
are integer and bigger than zero. The following example declares the first particle as a
j = 3/2 fermion and the second one as a d-boson.
```text
ParticleTypeSymmetries 4 5
```

### `ParticleTypeMaxNums`
t is possible to restrict the maximum number of each particle with this parameter. The
code takes this restriction into account during the construction of the quantum numbers
of the basis states. If a −1 for some particle is given, then the maximum number for that
particle will be the total particle number of the basis states (see TotalParticleNumber).
If a (integer) number k ≥ 0 is given, then the allowed particle number is 0,..,k. The
example gives no boundary for the first and third particle, but there can be only 0,1,2 or
3 particles of the second type.
```text
ParticleTypeMaxNums -1 3 -1
```

### `ParticleTypeSeniorities`
As ParticleTypeMaxNums but the seniorites given here will fix the seniorities to exactly
one value. This is useful for example for seniority conserving hamiltonians. Again a −1
means that there is no restriction. The seniorities are, of course, integer numbers bigger
than or equal zero.
**Example:**
```text
ParticleTypeSeniorities -1 -1 0
```

### `BasisLimit`
Using this parameter it is possible to give more specific restrictions to the basis, with
this parameter for example one can restrict the sum of the particle numbers of two or
more particle types to some integer number. Suppose the model has five distinguishable
particles. Then there have to be five integers number as arguments followed by either a
“=” or a “<=” followed again by an integer.
**Example:**
```text
Basislimit 0 1 0 2 1 <= 10
```
This example would restrict the basis in the following way: 0n_1 + 1n_2 + 0n_3 + 2n_4 + 1n_5 <= 10. where n_i is the particle number of the i-th particle type. For example this is useful to restrict the number of negative parity bosons for a
spdf-IBM1 calculation. There can be as many of such restrictions as needed.

### `TotalParticleNumber`
This would set the total particle number of the basis vectors to 12. Exactly one particle
number must be given. It is possible to do several calculations with one input file. The
example calculates all 0+ and all 2+ states but the total particle number will be fixed for
all calculations.
**Example:**
```text
TotalParticleNumber 12
```

### `TotalSpin`
Here the user can give a list of angular momenta in units of hbar/2. The hamiltonian matrix
will be constructed for each of these angular momenta 
**Example:**
```text
TotalSpin 0 4
```
This would calculate states with angular momenta 0 and 2.

### `TotalParity`
As TotalSpin but for the parity of each calculation. All parities must be either 1 or −1
and the number of parities has to be the same as the number of given angular momenta
with the parameter TotalSpin.
**Example:**
```text
TotalParity 1 -1
```

### `NumberOfStatesToCalc`
See also TotalSpin and TotalParity. For each calculation, i.e. for each given total spin and
parity, a number of states has to be given. The code will calculate for each Jπ always
the lowest number of states which is specified here. All states can be calculated with a
−1. The number of arguments to this parameter has to be the same as for TotalSpin and
TotalParity.
**Example:**
```text
NumberOfStatesToCalc -1 5
```

### `UseHamiltonian`
See also TotalSpin, TotalParity and NumberOfStatesToCalc. The code expects for each
calculation one string as arguments to this parameter. These strings can be choosen freely
and are just names of hamiltonians. These names can be the same if only one hamiltonian
is used. For each name the user has to give one hamiltonian explicitly in the formalism of
second quantization. For a more detailed explanation see H.
**Example:**
```text
UseHamiltonian H1 H2
```

### `WriteHamiltToFile`
This parameter is needed to perform a normal calculation. If it is given the arguments are
filenames (one for each calculation). The hamilton matrix will be written to these files.
The format will be the old Yale Sparse matrix format in ascii (see chapter 2.10 for an
introduction)
Example: WriteHamiltToFile h1.txt h2.txt 
**Example:**
```text
WriteHamiltToFile h1.txt h2.txt
```

### `DiagMethod`
Different algorithms for diagonalization are available. It depends on several properties of
the hamiltonian matrix which one is suited best. The user has to choose one for each
calculation by giving either a 1 or a 2. (A.L. Conley Note: program does not find method 2)
**Example:**
```text
DiagMethod 2 1
```

### `PrintBasis`
For each calculation it is possible to specify if the corresponding basis will be printed
out or not. This is done by using arguments ”y” for ”yes” and ”n” for ”no”. See also
PrintEigenvalues and PrintEigenvectors.
**Example:**
```text
PrintBasis y n
```

### `PrintEigenvalues`
The same as PrintBasis but for the eigenvalues. See also PrintEigenvectors.
**Example:**
```text
PrintEigenvalues y y
```

### `PrintEigenvectors`
Like PrintBasis and PrintEigenvalues but for the eigenvectors.
**Example:**
```text
PrintEigenvectors n n
```

### `PrintInputsummary`
If the argument of this parameter is ”y”, then an input summary will be printed out. The
user has to give a ”y” for ”yes” or a ”n” for ”no”.
**Example:**
```text
PrintInputsummary y
```

### `VerboseLevel`
There has to be one integer argument which must be 0, 1 or 2. If 0 is given, there will
be a low number of additional informations given during the calculations. The verbosity
increases for higher numbers. The given informations are helpful for real experts only.
**Example:**
```text
VerboseLevel 0
```

### `O` (Operator definition)
Definitions of for example transition operators are done with the parameter “O”. The
operator may have as many terms as needed but they all must be one body terms. The
first argument is an arbitrary string which is just a name. When specifying matrix ele-
ments which should be calculated one has to refer these names (see also M). The Example
defines an operator T2 := 0.4 ·(s† ×  ̃d)^(2) + 0.4 ·(d† × ̃s)^(2). Again, all angular momenta are in units of hbar/2.
**Example (defines T2):**
```text
O E2 s + d ~ 4 0.4 d + s ~ 4 0.4
```

### `M` (Matrix element request)
A line beginning with an ”M” is used to specify a reduced matrix element of a previously
defined operator which should be calculated (see also O). There can be as many of such
lines in the input file as needed. The example tells the program to calculate two reduced
matrix elements. The first one is the reduced matrix element of the operator with name
”T2” with respect to the third state of the first calculation and the second state of the
fourth calculation. The second matrix element is of the same operator but for the second
state of the second calculation and the second state of the fourth calculation.
**Example:**
```text
M 1 3 T2 4 2
M 2 2 T2 4 2
```

### `H` (Hamiltonian term)
Each line which starts with an ”H” gives a one or a two body term of a hamiltonian.
All of these terms will be added to get the final hamilton operator(s). The first argument
must be the name of the hamiltonian (see also UseHamiltonian). The user has to ensure
by himself that the hamiltonian is hermitian. The number of arguments differs for one
and two body terms.

- **One-body term**  
  **Example:**  
  The example gives a one body term 0.5 ·(d† ×  ̃d)(0) for the hamiltonian with the name “H3”
  ```text
  H H3 d + d ~ 0.5
  ```

- **Two-body terms**  
    The first example gives a two body term
    −2.2 ·[(s† ×  ̃d)^(2) ×(s† ×  ̃d)^(2)]^(0) for the hamiltonian with the name “name1” and
    the second one 0.3 ·[(s† × ̃s)^(2) ×(d† ×  ̃d)^(2)]^(0) for “H1”. If the names are the same,
    the terms will be added.
  **Example 1:**  
  ```text
  H name1 s + d ~ s + d ~ 4 -2.2
  ```

  **Example 2:**  
  
  ```text
  H H1 s + s + d ~ d ~ 0 0.3
  ```

---

**End of manual.**
