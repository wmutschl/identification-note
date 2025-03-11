clc; clear vars; clear globalvars; clearvars; close all;
dynare AnSchoModTheBuilder console -DINDEXATION=0 -DPREFSHOCK=0 -DMONPOL=1 -DMEASERR=1 -DVAROBSCOMBINATIONS=3 -DCOMPUTATIONS=1