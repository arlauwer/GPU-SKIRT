#!/bin/bash

# CMake
CMAKEPATH="$(which cmake)"

# Generate the build files
$CMAKEPATH -E make_directory ../debug
$CMAKEPATH -E chdir ../debug $CMAKEPATH \
	-DCMAKE_EXPORT_COMPILE_COMMANDS=1 \
	-DCMAKE_BUILD_TYPE=Debug \
	-DCMAKE_CXX_FLAGS="-g -O0" \
	-DCMAKE_EXE_LINKER_FLAGS="-rdynamic" \
	-L ../git
echo

# Perform the build
make -j ${1:-1} -C ../debug
echo
