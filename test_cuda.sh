#!/bin/bash
# Test if Cuda can be used for compiling

current_dir=`pwd`
script_dir=`dirname $0`
CUDA_VERSION=`nvcc --version 2>/dev/null | grep release | awk '{ print $(NF-1) }' | sed "s/,//g"`
if [[ $CUDA_VERSION != "" ]]; then
	printf "Using Cuda %s\n" $CUDA_VERSION >&2
else
	if [[ $DEVICE == "CUDA" ]]; then
		printf "Error: nvcc command does not exist/is not working properly.\n" >&2
	fi
	exit 1
fi
TARGETS_SUPPORTED=`nvcc --list-gpu-arch 2>/dev/null`
if [[ $TARGETS_SUPPORTED == "" ]]; then # might be an older Cuda version that doesn't have the --list-gpu-arch option
	TARGETS_SUPPORTED=`nvcc --help | grep -oP "compute_\d+" | sort -u`
fi
if [[ "$4" != "" ]]; then
	for T in $4; do
		SUPPORTED=`grep -o $T <<< $TARGETS_SUPPORTED`
		if [[ $SUPPORTED == "" ]]; then
			printf "Error: Specified compute target <$T> not supported by installed Cuda version.\n" >&2
			exit 1
		fi
	done
	TARGETS="$4"
else
	TARGETS=`awk -F'_' "{ if(\\$2>=$5) print \\$2 }" <<< "$TARGETS_SUPPORTED" | tr "\n" " "`
fi
printf "Compiling for targets: %s\n" "$TARGETS" >&2
cd "$script_dir"
if [[ ! -f "test_cuda" ]]; then
	$1 -I$2 -L$3 -lcuda -lcudart -o test_cuda test_cuda.cpp &> /dev/null
	test -e test_cuda && echo $TARGETS
else
	test -e test_cuda && echo $TARGETS
fi
cd "$current_dir"
