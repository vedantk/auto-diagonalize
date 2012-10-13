#!/bin/sh

if [ -z $1 ]; then
	echo "Usage: ./run.sh <source-file> [<flags>]"
	exit 1
fi

clang++ -emit-llvm $1 -o $1.bc -c
opt -load ./diagonalize.so -auto-diagonalize $1.bc -o $1.bc.diag
clang++ $1.bc.diag -o $1.out
