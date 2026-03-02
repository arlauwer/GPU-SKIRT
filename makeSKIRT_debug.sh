#!/bin/bash

# CMake
CMAKEPATH="$(which cmake)"

# Generate the build files
$CMAKEPATH -E make_directory ../debug
$CMAKEPATH -E chdir ../debug $CMAKEPATH -DCMAKE_EXPORT_COMPILE_COMMANDS=1 -DCMAKE_BUILD_TYPE:STRING=Debug -L ../git
echo

# Perform the build
make -j ${1:-1} -C ../debug
echo
