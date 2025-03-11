clc; clear vars; clear globalvars; clearvars; close all;
dynare AnSchoModTheBuilder console -DINDEXATION=1 -DPREFSHOCK=0 -DMONPOL=2 -DMEASERR=1 -DVAROBSCOMBINATIONS=3 -DCOMPUTATIONS=1