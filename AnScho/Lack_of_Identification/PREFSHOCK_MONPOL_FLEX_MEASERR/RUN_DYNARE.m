clc; clear vars; clear globalvars; clearvars; close all;
dynare AnSchoModTheBuilder console -DINDEXATION=0 -DPREFSHOCK=1 -DMONPOL=0 -DMEASERR=1 -DVAROBSCOMBINATIONS=3 -DCOMPUTATIONS=1