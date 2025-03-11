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
    jmodel = 1:10;
else
    allmodels = 0;
    jmodel = str2double(setting_models);        
end
    
%% Baseline SETTINGS
OPT.VAROBSCOMBINATIONS  = 2;
OPT.MSADJ               = 1;  % Always consider multisectoral adjustment costs
OPT.UTIL                = 0;
OPT.HABIT               = 0;
OPT.CAPUTIL             = 0;
OPT.LABOR               = 0;
OPT.INVESTSHOCK         = 0;
OPT.MONPOL              = 0;
OPT.COMPUTATIONS        = 1;  % Run Only Rank Checks of Identification

%% Case 1: BASELINE
if allmodels == 1 || jmodel == 1
    OPT1 = OPT;
    OPT1.namefolder = 'BASELINE';    
end

%% Case 2: CAPITAL UTILIZATION
if allmodels == 1 || jmodel == 2
    OPT2 = OPT;
    OPT2.CAPUTIL = 1;
    OPT2.namefolder = 'CAPUTILIZATION';
end

%% Case 3: INVESTMENT SPECIFIC TECHNOLOGICAL SHOCK
if allmodels == 1 || jmodel == 3
    OPT3 = OPT;
    OPT3.INVESTSHOCK = 1;
    OPT3.VAROBSCOMBINATIONS = 3;
    OPT3.namefolder = 'INVESTSPECSHOCK';
end

%% Case 4: EXTERNAL HABIT
if allmodels == 1 || jmodel == 4
    OPT4 = OPT;    
    OPT4.HABIT = 1;
    OPT4.namefolder = 'EXTERNALHABIT';
end

%% Case 5: INTERNAL HABIT
if allmodels == 1 || jmodel == 5
    OPT5 = OPT;
    OPT5.HABIT = 2;
    OPT5.namefolder = 'INTERNALHABIT';
end

%% Case 6: CRRA UTILITY NO HABIT
if allmodels == 1 || jmodel == 6
    OPT6 = OPT;
    OPT6.UTIL = 1;
    OPT6.namefolder = 'CRRA_NOHABIT';
end

%% Case 7: CRRA UTILITY EXTERNAL HABIT
if allmodels == 1 || jmodel == 7
    OPT7 = OPT;
    OPT7.UTIL = 1;
    OPT7.HABIT = 1;
    OPT7.namefolder = 'CRRA_EXTERNALHABIT';
end

%% Case 8: CRRA UTILITY INTERNAL HABIT
if allmodels == 1 || jmodel == 7
    OPT8 = OPT;
    OPT8.UTIL = 1;
    OPT8.HABIT = 2;
    OPT8.namefolder = 'CRRA_INTERNALHABIT';
end

%% Case 9: LABOR
if allmodels == 1 || jmodel == 9
    OPT9 = OPT;
    OPT9.LABOR = 1;
    OPT9.namefolder = 'LABOR';
end

%% Case 10: MONETARY POLICY
if allmodels == 1 || jmodel == 10
    OPT10 = OPT;
    OPT10.MONPOL = 1;
    OPT10.VAROBSCOMBINATIONS = 3;
    OPT10.namefolder = 'MONPOL';
end

%% Create Model Calls where we loop over iac=0,1,2 (intertemporal investment adjustment costs)
callstring = 'dynare KimModTheBuilder console -DIAC=%d -DMSADJ=%d -DUTIL=%d -DHABIT=%d -DCAPUTIL=%d -DINVESTSHOCK=%d -DLABOR=%d -DMONPOL=%d -DVAROBSCOMBINATIONS=%d -DCOMPUTATIONS=%d ';
for j=jmodel
    fprintf('Case %d\n',j);
    for jiac=1:2
        eval(sprintf('opt = OPT%d;',j));
        opt.IAC = jiac;
        if jiac == 1
            modelfolder = [opt.namefolder, '_IAC_LEVEL'];
        elseif jiac == 2
            modelfolder = [opt.namefolder, '_IAC_GROWTH'];
        end
        if exist(modelfolder, 'dir') == 7
            rmdir(modelfolder,'s');
            pause(1);
        end
        mkdir(modelfolder);
        pause(1);
        copyfile('../../utils/KimModTheBuilder.mod',modelfolder);
        copyfile('../../utils/LatexTable.m',modelfolder);
        calldynare = sprintf(callstring, opt.IAC, opt.MSADJ, opt.UTIL, opt.HABIT, opt.CAPUTIL, opt.INVESTSHOCK, opt.LABOR, opt.MONPOL, opt.VAROBSCOMBINATIONS, opt.COMPUTATIONS);
        fileID = fopen([modelfolder,'/RUN_DYNARE.m'],'w');
        fprintf(fileID,'clc; clear vars; clear globalvars; clearvars; close all;\n');
        fprintf(fileID,'%s',calldynare);
        fclose(fileID);
    end
end


%% Run Dynare for a model variant
fprintf('***********************************************************************************************************\n');
fprintf('Now you have to go into all model folders and run ''RUN_DYNARE.m'' with at least Dynare Version 4.6-unstable.\n');
fprintf('***********************************************************************************************************\n');