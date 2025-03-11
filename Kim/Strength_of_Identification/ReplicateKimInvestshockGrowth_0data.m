%% CLEAR WORKSPACE AND INITIALIZE
clc; clear vars; clear globalvars; clearvars; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
WEAKMODEL = 'InvestshockGrowth';
SIMULDATA = 1;
FIRSTOBS = 1234;
SAMPLESIZE = 50000;
HESSIAN = 0;
MCMC = 0;
MH_REPLIC = 1000000;
MH_NBLOCKS = 4;
CONSOLE = 1;
PARALLEL = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
RunKimStrength;
