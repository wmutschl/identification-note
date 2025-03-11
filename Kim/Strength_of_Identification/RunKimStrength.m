CLUSTER=0; %set to 1 for different parallel configuration file for cluster

mkdir([WEAKMODEL,'/',num2str(SAMPLESIZE)]);
oldfolder = cd([WEAKMODEL,'/',num2str(SAMPLESIZE)]);
copyfile('../../../../utils/KimModTheBuilder.mod');
calldynare = 'dynare KimModTheBuilder -DMSADJ=1 -DUTIL=0 -DHABIT=0 -DLABOR=0 -DMONPOL=0';
calldynare = sprintf('%s -DHESSIAN=%d -DMCMC=%d -DMH_REPLIC=%d -DMH_NBLOCKS=%d',calldynare,HESSIAN,MCMC,MH_REPLIC,MH_NBLOCKS);

switch WEAKMODEL
    case 'BaselineLevel'
        VAROBS = 'c';
        calldynare = [calldynare, ' -DIAC=1 -DCAPUTIL=0 -DINVESTSHOCK=0 -DVAROBS="' VAROBS, '"'];
        switch SAMPLESIZE
            case 100
                FINDMODEADVANCED = 0; OPTIMIZER = 4; JSCALEINTEGER = 0; JSCALEDECIMAL = 8;
            case 300
                FINDMODEADVANCED = 0; OPTIMIZER = 4; JSCALEINTEGER = 0; JSCALEDECIMAL = 7;
            case 900
                FINDMODEADVANCED = 0; OPTIMIZER = 4; JSCALEINTEGER = 0; JSCALEDECIMAL = 7;
            case 2700
                FINDMODEADVANCED = 0; OPTIMIZER = 4; JSCALEINTEGER = 0; JSCALEDECIMAL = 6;
            case 8100
                FINDMODEADVANCED = 0; OPTIMIZER = 4; JSCALEINTEGER = 0; JSCALEDECIMAL = 5;
            case 50000
                FINDMODEADVANCED = 0; OPTIMIZER = 4; JSCALEINTEGER = 0; JSCALEDECIMAL = 5;
        end
    case 'BaselineGrowth'
        VAROBS = 'c';
        calldynare = [calldynare, ' -DIAC=2 -DCAPUTIL=0 -DINVESTSHOCK=0 -DVAROBS="' VAROBS, '"'];
        switch SAMPLESIZE
            case 100
                FINDMODEADVANCED = 0; OPTIMIZER = 4; JSCALEINTEGER = 0; JSCALEDECIMAL = 8;
            case 300
                FINDMODEADVANCED = 0; OPTIMIZER = 4; JSCALEINTEGER = 0; JSCALEDECIMAL = 7;
            case 900
                FINDMODEADVANCED = 0; OPTIMIZER = 4; JSCALEINTEGER = 0; JSCALEDECIMAL = 7;
            case 2700
                FINDMODEADVANCED = 0; OPTIMIZER = 4; JSCALEINTEGER = 0; JSCALEDECIMAL = 7;
            case 8100
                FINDMODEADVANCED = 0; OPTIMIZER = 4; JSCALEINTEGER = 0; JSCALEDECIMAL = 6;
            case 50000
                FINDMODEADVANCED = 0; OPTIMIZER = 4; JSCALEINTEGER = 0; JSCALEDECIMAL = 6;
        end
    case 'InvestshockLevel'
        VAROBS = 'y c';
        calldynare = [calldynare, ' -DIAC=1 -DCAPUTIL=0 -DINVESTSHOCK=1 -DVAROBS="' VAROBS, '"'];
        switch SAMPLESIZE
            case 100
                FINDMODEADVANCED = 0; OPTIMIZER = 4; JSCALEINTEGER = 0; JSCALEDECIMAL = 8;
            case 300
                FINDMODEADVANCED = 0; OPTIMIZER = 4; JSCALEINTEGER = 0; JSCALEDECIMAL = 8;
            case 900
                FINDMODEADVANCED = 0; OPTIMIZER = 4; JSCALEINTEGER = 0; JSCALEDECIMAL = 7;
            case 2700
                FINDMODEADVANCED = 0; OPTIMIZER = 4; JSCALEINTEGER = 0; JSCALEDECIMAL = 7;
            case 8100
                FINDMODEADVANCED = 0; OPTIMIZER = 4; JSCALEINTEGER = 0; JSCALEDECIMAL = 7;
            case 50000
                FINDMODEADVANCED = 0; OPTIMIZER = 4; JSCALEINTEGER = 0; JSCALEDECIMAL = 7;
        end
    case 'InvestshockGrowth'
        VAROBS = 'y c';
        calldynare = [calldynare, ' -DIAC=2 -DCAPUTIL=0 -DINVESTSHOCK=1 -DVAROBS="' VAROBS, '"'];
        switch SAMPLESIZE
            case 100
                FINDMODEADVANCED = 0; OPTIMIZER = 4; JSCALEINTEGER = 0; JSCALEDECIMAL = 8;
            case 300
                FINDMODEADVANCED = 0; OPTIMIZER = 4; JSCALEINTEGER = 0; JSCALEDECIMAL = 8;
            case 900
                FINDMODEADVANCED = 0; OPTIMIZER = 4; JSCALEINTEGER = 0; JSCALEDECIMAL = 8;
            case 2700
                FINDMODEADVANCED = 0; OPTIMIZER = 4; JSCALEINTEGER = 0; JSCALEDECIMAL = 8;
            case 8100
                FINDMODEADVANCED = 0; OPTIMIZER = 4; JSCALEINTEGER = 0; JSCALEDECIMAL = 8;
            case 50000
                FINDMODEADVANCED = 0; OPTIMIZER = 4; JSCALEINTEGER = 0; JSCALEDECIMAL = 8;
        end
end
calldynare = sprintf([calldynare,' -DFINDMODEADVANCED=%d -DFIRSTOBS=%d -DSAMPLESIZE=%d -DOPTIMIZER=%d -DJSCALEINTEGER=%d -DJSCALEDECIMAL=%d'],FINDMODEADVANCED,FIRSTOBS,SAMPLESIZE,OPTIMIZER,JSCALEINTEGER,JSCALEDECIMAL);
if CONSOLE == 1
    calldynare = [calldynare, ' console'];
end
if PARALLEL == 1
    if CLUSTER == 1
        copyfile('../../../../utils/SetupForParallelCluster.ini');
        calldynare = [calldynare, ' conffile=''SetupForParallelCluster.ini'' parallel'];
    else
        copyfile('../../../../utils/SetupForParallel.ini');
        calldynare = [calldynare, ' conffile=''SetupForParallel.ini'' parallel'];
    end
end
% make sure workspace is clean
close all; clear globalvars;
clearvars -except oldfolder calldynare WEAKMODEL SAMPLESIZE SIMULDATA HESSIAN MCMC;
if SIMULDATA == 1
    eval([calldynare, ' -DCOMPUTATIONS=2']); % Simulate data
    clean_current_folder;                
    pause(3);
end
dlmwrite('currentmodelcall.txt',[calldynare, ' -DCOMPUTATIONS=3'],'delimiter','');
if HESSIAN == 1 || MCMC == 1
    eval([calldynare, ' -DCOMPUTATIONS=3']); % Run estimation
    if HESSIAN == 1
        save(['../weakresults_hessian_', num2str(SAMPLESIZE)],'weakresults_hessian');
    end
    if MCMC == 1
        save(['../weakresults_mcmc_', num2str(SAMPLESIZE)],'weakresults_mcmc');
    end
end
% Clean up
cd(oldfolder);
