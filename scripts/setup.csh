#!/bin/csh
#source /cvmfs/sft.cern.ch/lcg/views/LCG_89/x86_64-centos7-gcc7-opt/setup.csh
#source /cvmfs/sft.cern.ch/lcg/views/ROOT-latest/x86_64-centos7-gcc7-opt/setup.csh

setenv LD_LIBRARY_PATH ./lib:DynamicTTree/lib/:CfgManager/lib/:$LD_LIBRARY_PATH
setenv ROOT_INCLUDE_PATH ./interface:DynamicTTree/interface/:CfgManager/interface/:$ROOT_INCLUDE_PATH
