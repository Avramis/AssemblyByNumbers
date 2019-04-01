#!/bin/bash
echo '#############################################'
echo '#                                           #'
echo '#               Intel Compiler              #'
echo '# Compiling AssemblyByNumbers initiation... #'
echo '#                                           #'


cwd=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"/Code/"
cd $cwd


icpc -Wall -O3 -std=c++11 -o AssemblyByNumbers *.cpp *.c
mv AssemblyByNumbers ../Exec

echo '#                                           #'
echo '# Compiling AssemblyByNumbers completed...  #'
echo '#                                           #'
echo '#############################################'

