%% CLEAR WORKSPACE AND INITIALIZE FOLDERS
clc; clear vars; clear globalvars; clearvars; close all;

for strength_table_model = [{'Baseline'}, {'MonPolSteadyStateGap'}, {'Indexation'}, {'Prefshock'}]

    %% Create Strength of Identification Tables
    oldfolder = cd(strength_table_model{:});
    copyfile('../../../utils/LatexTable.m')
    mkdir('tables')
    SAMPLESIZE = [100; 300; 900; 2700; 8100];
    load('100/AnSchoModTheBuilder_results.mat','bayestopt_');
    % make table 1
    inputTbl.tableColLabels = bayestopt_.name';
    inputTbl.tableRowLabels = {'100' '300' '900' '2700' '8100'};
    inputTbl.booktabs = 1;
    inputTbl.longtable = 1;
    inputTbl.dataFormat = {'%.3f'};
    inputTbl.makeCompleteLatexDocument = 1;
    WEAKRESULTS.hessian = zeros(5,size(bayestopt_.name,1));
    WEAKRESULTS.mcmc    = zeros(5,size(bayestopt_.name,1));
    % make table 2
    inputTbl2.tableColLabels = bayestopt_.name';
    inputTbl2.tableRowLabels = {'300/100' '900/300' '2700/900' '8100/2700'};
    inputTbl2.booktabs = 1;
    inputTbl2.longtable = 1;
    inputTbl2.dataFormat = {'%.3f'};
    inputTbl2.makeCompleteLatexDocument = 1;
    WEAKRESULTS2.hessian = zeros(4,size(bayestopt_.name,1));
    WEAKRESULTS2.mcmc    = zeros(4,size(bayestopt_.name,1));    
    for jj = 1:5
        load(sprintf('weakresults_hessian_%d.mat',SAMPLESIZE(jj)),'weakresults_hessian');
        load(sprintf('weakresults_mcmc_%d.mat',SAMPLESIZE(jj)),'weakresults_mcmc');             
        WEAKRESULTS.hessian(jj,:) = weakresults_hessian;
        WEAKRESULTS.mcmc(jj,:) = weakresults_mcmc;
        if jj>1
            WEAKRESULTS2.hessian(jj-1,:) = ((WEAKRESULTS.hessian(jj-1,:).^(-1))./SAMPLESIZE(jj-1)) ./ ((WEAKRESULTS.hessian(jj,:).^(-1))./SAMPLESIZE(jj));
            WEAKRESULTS2.mcmc(jj-1,:) = ((WEAKRESULTS.mcmc(jj-1,:).^(-1))./SAMPLESIZE(jj-1)) ./ ((WEAKRESULTS.mcmc(jj,:).^(-1))./SAMPLESIZE(jj));
        end   
    end
    for strres = {'hessian', 'mcmc'}            
        inputTbl.data = WEAKRESULTS.(strres{:});
        inputTbl.tableCaption = ['Bayesian Weak Identification An Schorfheide ', strres{:} ' method'];
        inputTbl.tableLabel = ['tbl:WeakAnScho_', strres{:}];
        filenam = ['table_AnSchoWeak_', strres{:},'.tex'];
        latex = LatexTable(inputTbl);
        fid=fopen(['tables/' filenam],'w');
        [nrows,ncols] = size(latex);
        for row = 1:nrows
            fprintf(fid,'%s\n',latex{row,:});
        end
        fclose(fid);
        oldfolder2 = cd('tables');
        for kk = 1:3
            system(['pdflatex -halt-on-error ' filenam])
        end
        cd(oldfolder2);
    end
    
    for strres = {'hessian', 'mcmc'}            
        inputTbl2.data = WEAKRESULTS2.(strres{:});
        inputTbl2.tableCaption = ['Bayesian Weak Identification An Schorfheide Convergence Ratios', strres{:} ' method'];
        inputTbl2.tableLabel = ['tbl:WeakAnSchoConvergenceRatios_', strres{:}];
        filenam = ['table_AnSchoWeak_ConvRatios_', strres{:},'.tex'];
        latex = LatexTable(inputTbl2);
        fid=fopen(['tables/' filenam],'w');
        [nrows,ncols] = size(latex);
        for row = 1:nrows
            fprintf(fid,'%s\n',latex{row,:});
        end
        fclose(fid);
        oldfolder2 = cd('tables');
        for kk = 1:3
            system(['pdflatex -halt-on-error ' filenam])
        end
        cd(oldfolder2);
    end
    
    cd(oldfolder);
    fprintf('************************************************************************************************************\n');

end