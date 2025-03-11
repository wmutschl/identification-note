CLUSTER=0; %set to 1 for different parallel configuration file for cluster

mkdir([WEAKMODEL,'/',num2str(SAMPLESIZE)]);
oldfolder = cd([WEAKMODEL,'/',num2str(SAMPLESIZE)]);
copyfile('../../../../utils/AnSchoModTheBuilder.mod');
calldynare = sprintf('dynare AnSchoModTheBuilder -DHESSIAN=%d -DMCMC=%d -DMH_REPLIC=%d -DMH_NBLOCKS=%d',HESSIAN,MCMC,MH_REPLIC,MH_NBLOCKS);

switch WEAKMODEL
    case 'Baseline'
        VAROBS = 'YGR INFL INT';
        calldynare = [calldynare, ' -DMONPOL=0 -DINDEXATION=0 -DPREFSHOCK=0 -DMEASERR=0 -DVAROBS="' VAROBS, '"'];
        switch SAMPLESIZE
            case 100
                FINDMODEADVANCED = 1; OPTIMIZER = 4; JSCALEINTEGER = 0; JSCALEDECIMAL = 55;
            case 300
                FINDMODEADVANCED = 1; OPTIMIZER = 4; JSCALEINTEGER = 0; JSCALEDECIMAL = 6;
            case 900
                FINDMODEADVANCED = 1; OPTIMIZER = 4; JSCALEINTEGER = 0; JSCALEDECIMAL = 65;
            case 2700
                FINDMODEADVANCED = 1; OPTIMIZER = 4; JSCALEINTEGER = 0; JSCALEDECIMAL = 65;
            case 8100
                FINDMODEADVANCED = 1; OPTIMIZER = 4; JSCALEINTEGER = 0; JSCALEDECIMAL = 75;
            case 50000
                FINDMODEADVANCED = 1; OPTIMIZER = 4; JSCALEINTEGER = 0; JSCALEDECIMAL = 75;
        end
    case 'MonPolSteadyStateGap'
        VAROBS = 'YGR INFL INT';
        calldynare = [calldynare, ' -DMONPOL=1 -DINDEXATION=0 -DPREFSHOCK=0 -DMEASERR=0 -DVAROBS="' VAROBS, '"'];
        switch SAMPLESIZE
            case 100
                FINDMODEADVANCED = 0; OPTIMIZER = 4; JSCALEINTEGER = 0; JSCALEDECIMAL = 6;
            case 300
                FINDMODEADVANCED = 0; OPTIMIZER = 4; JSCALEINTEGER = 0; JSCALEDECIMAL = 6;
            case 900
                FINDMODEADVANCED = 0; OPTIMIZER = 4; JSCALEINTEGER = 0; JSCALEDECIMAL = 6;
            case 2700
                FINDMODEADVANCED = 0; OPTIMIZER = 4; JSCALEINTEGER = 0; JSCALEDECIMAL = 6;
            case 8100
                FINDMODEADVANCED = 0; OPTIMIZER = 4; JSCALEINTEGER = 0; JSCALEDECIMAL = 6;
            case 50000
                FINDMODEADVANCED = 0; OPTIMIZER = 4; JSCALEINTEGER = 0; JSCALEDECIMAL = 6;
        end
    case 'Indexation'
        VAROBS = 'YGR INFL INT';
        calldynare = [calldynare, ' -DMONPOL=0 -DINDEXATION=1 -DPREFSHOCK=0 -DMEASERR=0 -DVAROBS="' VAROBS, '"'];
        switch SAMPLESIZE
            case 100
                FINDMODEADVANCED = 1; OPTIMIZER = 4; JSCALEINTEGER = 0; JSCALEDECIMAL = 55;
            case 300
                FINDMODEADVANCED = 1; OPTIMIZER = 4; JSCALEINTEGER = 0; JSCALEDECIMAL = 55;
            case 900
                FINDMODEADVANCED = 1; OPTIMIZER = 4; JSCALEINTEGER = 0; JSCALEDECIMAL = 65;
            case 2700
                FINDMODEADVANCED = 1; OPTIMIZER = 4; JSCALEINTEGER = 0; JSCALEDECIMAL = 65;
            case 8100
                FINDMODEADVANCED = 1; OPTIMIZER = 4; JSCALEINTEGER = 0; JSCALEDECIMAL = 75;
            case 50000
                FINDMODEADVANCED = 1; OPTIMIZER = 4; JSCALEINTEGER = 0; JSCALEDECIMAL = 75;
        end
    case 'Prefshock'
        VAROBS = 'YGR INFL INT';
        calldynare = [calldynare, ' -DMONPOL=0 -DINDEXATION=0 -DPREFSHOCK=1 -DMEASERR=0 -DVAROBS="' VAROBS, '"'];
        switch SAMPLESIZE
            case 100
                FINDMODEADVANCED = 1; OPTIMIZER = 4; JSCALEINTEGER = 0; JSCALEDECIMAL = 4;
            case 300
                FINDMODEADVANCED = 1; OPTIMIZER = 4; JSCALEINTEGER = 0; JSCALEDECIMAL = 5;
           case 900
                FINDMODEADVANCED = 1; OPTIMIZER = 4; JSCALEINTEGER = 0; JSCALEDECIMAL = 6;
            case 2700
                FINDMODEADVANCED = 1; OPTIMIZER = 4; JSCALEINTEGER = 0; JSCALEDECIMAL = 5;
            case 8100
                FINDMODEADVANCED = 1; OPTIMIZER = 4; JSCALEINTEGER = 0; JSCALEDECIMAL = 8;
            case 50000
                FINDMODEADVANCED = 1; OPTIMIZER = 4; JSCALEINTEGER = 0; JSCALEDECIMAL = 8;
        end
end
calldynare = sprintf([calldynare,' -DFINDMODEADVANCED=%d -DFIRSTOBS=%d -DSAMPLESIZE=%d -DOPTIMIZER=%s -DJSCALEINTEGER=%d -DJSCALEDECIMAL=%d'],FINDMODEADVANCED,FIRSTOBS,SAMPLESIZE,mat2str(OPTIMIZER),JSCALEINTEGER,JSCALEDECIMAL);
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
