#!/bin/bash

# Enter directory in which the setup.sh file is located.
cd $(dirname $(readlink -f "$0")) &&

# Initialize all submodules to be up to date with GitHub versions.
git submodule update --init --recursive &&

# Delete old build directory and re-install the library using cmake.
rm -rf build && cmake -E make_directory build && cd build && cmake ..
