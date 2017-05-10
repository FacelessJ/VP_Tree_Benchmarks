#!/bin/bash

buildType="Release"
if [[ $# -ge 1 ]] ; then
	if [[ $1 == "Debug" ]]; then
		buildType="Debug"
	fi
fi

mkdir -p build_linux
cd build_linux
cmake -DCMAKE_BUILD_TYPE=$buildType .. && cmake --build .
cd ..
