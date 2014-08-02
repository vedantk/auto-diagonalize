#!/bin/bash

if [ -z $1 ]; then
	echo "Usage: ./run.sh <source-file> [<flags>]"
	exit 1
fi

clang++ $2 -std=c++11 -emit-llvm $1 -o $1.bc -c
opt -load ./diagonalize.so -auto-diagonalize $1.bc -o $1.bc.diag
llc $1.bc.diag -o $1.bc.diag.s
gcc $1.bc.diag.s -lm -o $1.out
