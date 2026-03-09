#!/bin/bash

# CMake
CMAKEPATH="$(which cmake)"

# Generate the build files
$CMAKEPATH -E make_directory ../release
$CMAKEPATH -E chdir ../release $CMAKEPATH \
	-DCMAKE_EXPORT_COMPILE_COMMANDS=1 \
	-DCMAKE_BUILD_TYPE=Release \
	-L ../git
echo

# Perform the build
make -j ${1:-1} -C ../release
echo
