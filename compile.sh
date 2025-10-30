#!/bin/sh
set -e  # stop immediately if any command fails

# --- Compiler and optimization flags ---
CXX="${CXX:-g++}"
CXXFLAGS="${CXXFLAGS:--O3 -std=c++17 -march=native -Wall -Wextra}"

# Eigen (header-only) lives here on Apple Silicon
INCLUDES="-I/opt/homebrew/include/eigen3"

# Link to Accelerate for optimized BLAS/LAPACK (optional but fast)
LIBS="-framework Accelerate -lpthread -ldl -lm"

cd ./Classes

echo "compiling csubroutines.cpp .."
rm -f ./csubroutines.o
$CXX $CXXFLAGS $INCLUDES -c ./csubroutines.cpp

echo "compiling cdrecouple.cpp .."
rm -f ./cdrecouple.o
$CXX $CXXFLAGS $INCLUDES -c ./cdrecouple.cpp

echo "compiling carpackwrapper.cpp .."
rm -f ./carpackwrapper.o
$CXX $CXXFLAGS $INCLUDES -c ./carpackwrapper.cpp

echo "compiling cidparticlestateqn2.cpp .."
rm -f ./cidparticlestateqn2.o
$CXX $CXXFLAGS $INCLUDES -c ./cidparticlestateqn2.cpp

echo "compiling clinearbasisqn2.cpp .."
rm -f ./clinearbasisqn2.o
$CXX $CXXFLAGS $INCLUDES -c ./clinearbasisqn2.cpp

echo "compiling cdreducedmatrix2.cpp .."
rm -f ./cdreducedmatrix2.o
$CXX $CXXFLAGS $INCLUDES -c ./cdreducedmatrix2.cpp

echo "compiling cdidparticlecfp.cpp .."
rm -f ./cdidparticlecfp.o
$CXX $CXXFLAGS $INCLUDES -c ./cdidparticlecfp.cpp

cd ..

echo "linking arbmodel .."
$CXX $CXXFLAGS -o arbmodel \
    ./arbmodel.c \
    ./Classes/cdreducedmatrix2.o \
    ./Classes/clinearbasisqn2.o \
    ./Classes/cidparticlestateqn2.o \
    ./Classes/cdidparticlecfp.o \
    ./Classes/cdrecouple.o \
    ./Classes/csubroutines.o \
    $LIBS -x c++

echo "✅ Build complete: ./arbmodel"