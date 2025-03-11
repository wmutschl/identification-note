% CLEAR WORKSPACE AND INITIALIZE FOLDERS
clc; clear vars; clear globalvars; clearvars; close all;

% Create Model Variants (if no option is set, then it creates all models new
fprintf('************************************************************************************************************\n');
fprintf('All the models have been run already, so you can just check the log files in the model folders.\n');
fprintf('The identification results are stored in a mat file called ''RESULTS_TABLE.mat''.\n');
fprintf('Latex Tables are also generated and inside the model folders.\n');

setting_models = input('To recreate all model files enter ''all''\nTo recreate a specific model variant enter the number as string e.g. ''1'': ');
if strcmp(setting_models,'all')
    allmodels = 1;
    jmodel = 1:4;
else
    allmodels = 0;
    jmodel = str2double(setting_models);
end

%% BASELINE SETTINGS
OPT.VAROBSCOMBINATIONS  = 3;
OPT.INDEXATION          = 0;
OPT.PREFSHOCK           = 0;
OPT.COMPUTATIONS        = 1;

%% Case 1: BASELINE
if allmodels == 1 || jmodel == 1
    OPT1 = OPT;
    OPT1.namefolder = 'BASELINE'; 
end

%% Case 2: INDEXATION
if allmodels == 1 || jmodel == 2
    OPT2 = OPT;
    OPT2.INDEXATION = 1;
    OPT2.namefolder = 'INDEXATION'; 
end

%% Case 3: PREFSHOCK
if allmodels == 1 || jmodel == 3
    OPT3 = OPT;
    OPT3.PREFSHOCK = 1;
    OPT3.namefolder = 'PREFSHOCK'; 
end

%% Case 4: INDEXATION AND PREFSHOCK
if allmodels == 1 || jmodel == 4
    OPT4 = OPT;
    OPT4.INDEXATION = 1;
    OPT4.PREFSHOCK = 1;
    OPT4.namefolder = 'INDEXATION_AND_PREFSHOCK'; 
end


%% Create Model Calls, we loop over jmp=0,1,2 (monetary policy rule) and jmeaserr=0,1 (measurement errors)
callstring = 'dynare AnSchoModTheBuilder console -DINDEXATION=%d -DPREFSHOCK=%d -DMONPOL=%d -DMEASERR=%d -DVAROBSCOMBINATIONS=%d -DCOMPUTATIONS=%d';

for j=jmodel
    fprintf('Case %d\n',j);
    for jmeaserr = 0:1
        for jmp=0:3
            eval(sprintf('opt = OPT%d;',j));
            opt.MEASERR = jmeaserr;
            opt.MONPOL = jmp;
            if jmp == 0
                modelfolder = [opt.namefolder, '_MONPOL_FLEX'];
            elseif jmp == 1
                modelfolder = [opt.namefolder, '_MONPOL_STEADYSTATE'];
            elseif jmp == 2
                modelfolder = [opt.namefolder, '_MONPOL_GROWTH'];
            elseif jmp == 3
                modelfolder = [opt.namefolder, '_MONPOL_SW'];
            end
            if jmeaserr == 1
                modelfolder = [modelfolder, '_MEASERR'];
            end
            if exist(modelfolder, 'dir') == 7
                rmdir(modelfolder,'s');
                pause(1);
            end
            mkdir(modelfolder);
            pause(1);
            copyfile('../../utils/AnSchoModTheBuilder.mod',modelfolder);
            copyfile('../../utils/LatexTable.m',modelfolder);
            calldynare = sprintf(callstring, opt.INDEXATION, opt.PREFSHOCK, opt.MONPOL, opt.MEASERR, opt.VAROBSCOMBINATIONS, opt.COMPUTATIONS);
            fileID = fopen([modelfolder,'/RUN_DYNARE.m'],'w');
            fprintf(fileID,'clc; clear vars; clear globalvars; clearvars; close all;\n');
            fprintf(fileID,'%s',calldynare);
            fclose(fileID);
        end
    end
end

%% Run Dynare for a model variant
fprintf('***********************************************************************************************************\n');
fprintf('Now you have to go into all model folders and run ''RUN_DYNARE.m'' with at least Dynare Version 4.6-unstable.\n');
fprintf('***********************************************************************************************************\n');