#!/bin/csh
source /cvmfs/sft.cern.ch/lcg/views/LCG_89/x86_64-slc6-gcc62-opt/setup.csh
source /cvmfs/sft.cern.ch/lcg/releases/ROOT/6.10.02-19565/x86_64-slc6-gcc62-opt/bin/thisroot.csh

setenv LD_LIBRARY_PATH ./lib:DynamicTTree/lib/:CfgManager/lib/:$LD_LIBRARY_PATH
setenv ROOT_INCLUDE_PATH ./interface:DynamicTTree/interface/:CfgManager/interface/:$ROOT_INCLUDE_PATH
