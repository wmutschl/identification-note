%% CLEAR WORKSPACE AND INITIALIZE
clc; clear vars; clear globalvars; clearvars; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
WEAKMODEL = 'InvestshockGrowth';
SIMULDATA = 0;
FIRSTOBS = 1234;
SAMPLESIZE = 100;
HESSIAN = 1;
MCMC = 1;
MH_REPLIC = 1000000;
MH_NBLOCKS = 4;
CONSOLE = 1;
PARALLEL = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
RunKimStrength;
