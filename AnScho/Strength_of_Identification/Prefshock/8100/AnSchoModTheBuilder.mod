% Replication files for:
% Ivashchenko and Mutschler - On the effect of observables, functional specifications, 
% and modeling choice on local identification in linearized DSGE models
% =========================================================================
% Note that this Mod file needs to be run with at least Dynare 4.6-unstable
% =========================================================================
% This mod file implements
% (1) several different features into the model of An and Schorfheide (2007) - Bayesian Analysis of DSGE Models; namely:
%     - Preference shock in the shape of a discount factor shifter: with or without
%     - Partial indexation scheme: full or partial indexation as in Smets and Wouters (2007)
%     - Different monetary policy rules: different definition of output gap:
%       deviations from (i) flex-price output, (ii) steady state output or (iii) output-growth
%
% (2) lack of identification tests of Iskrev (2010), Komunjer and Ng (2011)
%     and Qu and Tkachenko (2012) for all possible subsets of observations.
%     The intermediate results are saved into structures and a latex table
%     is compiled.
%
% (3) strength of identification test (Bayesian Learning Rate indicator)
%     of Koop, Pesaran, Smith (2013). The average posterior precision
%     for growing subsamples is computed based on both the Hessian as well 
%     as a full-fledged MCMC.
% -------------------------------------------------------------------------
% Copyright (C) 2019 Willi Mutschler (willi@mutschler.eu)
% 
% This is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% It is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% For a copy of the GNU General Public License,
% see <http://www.gnu.org/licenses/>.
% =========================================================================

% =========================================================================
% Settings For Dynare Preprocessor
% =========================================================================
% - PREFSHOCK:          0: Model without a discount factor shifter
%                       1: Model with a discount factor shifter
% -------------------------------------------------------------------------
% - INDEXATION:         0: Model with full indexation
%                       1: Model with partial indexation
% -------------------------------------------------------------------------
% - MONPOL:             Monetary policy rule that reacts to deviations of
%                       0: output from flex-price output
%                       1: output from steady state output
%                       2: output growth from steady state output growth
%                       3: output from flex-price output and output growth from flex-price output growth
% -------------------------------------------------------------------------
% - MEASERR:            0: Model without measurement errors
%                       1: Model with measurement errors on YGR, INFL and INT
% -------------------------------------------------------------------------
% - COMPUTATIONS:       0: Compute steady state, check Blanchard-Khan conditions, do model diagnostics
%                       1: Lack of Identification Analysis
%                       2: Simulate data
%                       3: Strength of Identification Analysis
% -------------------------------------------------------------------------
% - VAROBSCOMBINATIONS: [integer] maximum dimension of all possible subsets
%                       of observable variables to consider (only if COMPUTATIONS == 1)
% -------------------------------------------------------------------------
% - SAMPLESIZE:         [integer] sample size for data simulation (only if COMPUTATIONS >= 2)
% -------------------------------------------------------------------------
% - VAROBS:             [string] observable variables (only if COMPUTATIONS == 3)
% -------------------------------------------------------------------------
% - FIRSTOBS:           [integer] selects first observation (only if COMPUTATIONS == 3)
% -------------------------------------------------------------------------
% - HESSIAN:            1: compute mode and hessian at the mode (only if COMPUTATIONS == 3)
%                       0: skip mode computations (only if COMPUTATIONS == 3)
% -------------------------------------------------------------------------
% - FINDMODEADVANCED:   0: Use only one optimizer for mode finding (only if COMPUTATIONS == 3)
%                       1: Advanced mode finding procedure which loops over different
%                          optimizers and then uses Dynare's Monte Carlo "Optimizer" for a
%                          well-behaved Hessian in the relevant paramter
%                          space (only if COMPUTATIONS == 3)
% -------------------------------------------------------------------------
% - OPTIMIZER:          [integer] select which optimizer to use for mode finding
%                       (only if COMPUTATIONS == 3 && FINDMODEADVANCED == 0)
% -------------------------------------------------------------------------
% - MCMC:               1: compute posterior using metropolis hastings algorithm
%                       0: skip MCMC
%                       (only if COMPUTATIONS == 3)
% -------------------------------------------------------------------------
% - MH_NBLOCKS:         [integer] number of MCMC chains (only if COMPUTATIONS == 3)
% -------------------------------------------------------------------------
% - MH_REPLIC:          [integer] number of MCMC draws in each chain (only if COMPUTATIONS == 3)
% -------------------------------------------------------------------------
% - JSCALEINTEGER:      [integer] integer part of jscale tuning parameter (e.g. the 1 in 1.7) (only if COMPUTATIONS == 3)
% -------------------------------------------------------------------------
% - JSCALEDECIMAL:      [integer] decimal part of jscale tuning parameter (e.g. the 7 in 1.7) (only if COMPUTATIONS == 3)
% -------------------------------------------------------------------------

% =========================================================================
% Declare endogenous variables
% =========================================================================
var
YGR         ${YGR}$             (long_name='output growth rate (quarter-on-quarter)')
INFL        ${INFL}$            (long_name='annualized inflation rate')
INT         ${INT}$             (long_name='annualized nominal interest rate')
y           ${y}$               (long_name='detrended output (Y/A)')
c           ${c}$               (long_name='detrended consumption (C/A)')
r           ${R}$               (long_name='nominal interest rate')
p           ${\pi}$             (long_name='gross inflation rate')
g           ${g}$               (long_name='government consumption process (g = 1/(1-G/Y))')
z           ${z}$               (long_name='shifter to steady-state technology growth')
@#if PREFSHOCK == 1
zeta        ${\zeta}$           (long_name='shifter to discount factor')
@#endif
; % [var] end


% =========================================================================
% Declare exogenous variables
% =========================================================================
varexo
epsr        ${\varepsilon^R}$       (long_name='monetary policy shock')
epsg        ${\varepsilon^g}$       (long_name='government spending shock')
epsz        ${\varepsilon^z}$       (long_name='total factor productivity growth shock')
@#if PREFSHOCK == 1
epszeta     ${\varepsilon^\zeta}$   (long_name='discount factor shifter shock')
@#endif
@#if MEASERR == 1
measerr_YGR     ${\epsilon^{YGR}}$   (long_name='measurement error on YGR')
measerr_INFL    ${\epsilon^{INFL}}$  (long_name='measurement error on INFL')
measerr_INT     ${\epsilon^{INT}}$   (long_name='measurement error on INT')
@#endif
; % [varexo] end


% =========================================================================
% Declare parameters
% =========================================================================
parameters
RA          ${r_{A}}$               (long_name='annualized steady-state real interest rate')
PA          ${\pi^{(A)}}$           (long_name='annualized target inflation rate')
GAMQ        ${\gamma^{(Q)}}$        (long_name='quarterly steady-state growth rate of technology')
TAU         ${\tau}$                (long_name='inverse of intertemporal elasticity of subsitution')
NU          ${\nu}$                 (long_name='inverse of elasticity of demand in Dixit Stiglitz aggregator')
PSIP        ${\psi_\pi}$            (long_name='Taylor rule sensitivity parameter to inflation deviations')
PSIY        ${\psi_y}$              (long_name='Taylor rule sensitivity parameter to output deviations')
@#if MONPOL == 3
PSIDY       ${\psi_{dy}}$           (long_name='Taylor rule sensitivity parameter to trend output and growth deviations')
@#endif
RHOR        ${\rho_R}$              (long_name='Taylor rule persistence')
RHOG        ${\rho_{g}}$            (long_name='persistence government spending process')
RHOZ        ${\rho_z}$              (long_name='persistence TFP growth rate process')
SIGR        ${\sigma_R}$            (long_name='standard deviation monetary policy shock')
SIGG        ${\sigma_{g}}$          (long_name='standard deviation government spending process')
SIGZ        ${\sigma_z}$            (long_name='standard deviation TFP growth shock')
PHI         ${\phi}$                (long_name='Rotemberg adjustment cost parameter')
C_o_Y       ${\bar{C}/\bar{Y}}$     (long_name='steady state consumption to output ratio')
@#if INDEXATION == 1
IOTAP       ${\iota_p}$             (long_name='partial inflation indexation')
@#endif
@#if PREFSHOCK == 1
RHOZETA     ${\rho_\zeta}$          (long_name='persistence discount rate shifter')
SIGZETA     ${\sigma_\zeta}$        (long_name='standard deviation discount rate shifter shock')
@#endif
@#if MEASERR == 1
SIGYGR      ${\sigma_{YGR}}$        (long_name='standard deviation measurement errors YGR')
SIGINFL     ${\sigma_{INFL}}$       (long_name='standard deviation measurement errors INFL')
SIGINT      ${\sigma_{INT}}$        (long_name='standard deviation measurement errors INT')
@#endif
; % parameter block end


% =========================================================================
% Calibrate parameter values
% =========================================================================
% Relevant calibration for full model as in An and Schorfheide (2007)
RA      = 1;
PA      = 3.2;
GAMQ    = 0.55;
TAU     = 2;
NU      = 0.1;
PSIP    = 1.5;
PSIY    = 0.125;
PSIDY   = 0.2;
RHOR    = 0.75;
RHOG    = 0.95;
RHOZ    = 0.9;
SIGR    = 0.2;
SIGG    = 0.6;
SIGZ    = 0.3;
C_o_Y   = 0.85;
PHI     = 50;
@#if INDEXATION == 1
IOTAP   = 0.5;
@#endif
@#if PREFSHOCK == 1
RHOZETA = 0.75;
SIGZETA = 0.2;
@#endif
@#if MEASERR == 1
SIGYGR  = 0.23;
SIGINFL = 0.56;
SIGINT  = 0.66;
@#endif


% =========================================================================
% Model equations
% =========================================================================
model;
% -------------------------------------------------------------------------
% Auxiliary parameters and variables
#GAM = 1+GAMQ/100;
#BET = 1/(1+RA/400);
#PSTAR = 1+PA/400;
#G = 1/C_o_Y;
%Flex-price output
#ystar = (1-NU)^(1/TAU)*g;
#ystarback = (1-NU)^(1/TAU)*g(-1);
% -------------------------------------------------------------------------
% Monetary policy specification
@#if MONPOL == 0
% Flex-price output gap definition
#RSTAR = steady_state(r)*(p/PSTAR)^PSIP*(y/ystar)^PSIY;
@#endif
@#if MONPOL == 1
% Steady-state output gap definition
#RSTAR = steady_state(r)*(p/PSTAR)^PSIP*(y/steady_state(y))^PSIY;
@#endif
@#if MONPOL == 2
% Output growth definition
#RSTAR = steady_state(r)*(p/PSTAR)^PSIP*(y/y(-1)*z)^PSIY;
@#endif
@#if MONPOL == 3
% Trend output and growth definition
#RSTAR = steady_state(r)*(p/PSTAR)^PSIP*(y/ystar)^PSIY*( (y/y(-1))/(ystar/ystarback) )^PSIDY;
@#endif
% -------------------------------------------------------------------------
% Indexation rule
@#if INDEXATION == 1
#gammap = p^IOTAP*steady_state(p)^(1-IOTAP);
#gammapback = p(-1)^IOTAP*steady_state(p)^(1-IOTAP);
@#else
#gammap = steady_state(p);
#gammapback = steady_state(p);
@#endif
% -------------------------------------------------------------------------
% Auxiliary expression for preference shock
@#if PREFSHOCK == 0
#zeta = 1;
#zetap = 1;
@#else
#zetap = zeta(+1);
@#endif
% -------------------------------------------------------------------------
% Marginal utility and foc wrt c
#du_dc = c^(-TAU);
#dup_dcp = c(+1)^(-TAU);
#lam = du_dc*zeta;
#lamp = dup_dcp*zetap;

% -------------------------------------------------------------------------
% Actual Model Equations
% -------------------------------------------------------------------------
[name='Euler equation']
lam = BET*lamp*r/p(+1)/(GAM*z(+1));

[name='Price setting based on Rotemberg quadratic price adjustment costs and Dixit/Stiglitz aggregator']
1 = 1/NU*(1-(lam/zeta)^(-1)) + PHI*(p-gammapback)*p - PHI/(2*NU)*(p-gammapback)^2 - PHI*BET*(lamp/lam * y(+1)/y * (p(+1)-gammap)*p(+1));

[name='Market clearing: aggregate supply equals aggregate demand']
y - PHI/2*(p-gammapback)^2*y = c + (1-1/g)*y;

[name='Taylor rule']
r = RSTAR^(1-RHOR)*r(-1)^RHOR*exp(SIGR/100*epsr);

[name='Government spending process']
log(g) = (1-RHOG)*log(G) + RHOG*log(g(-1)) + SIGG/100*epsg;

[name='Technology growth process']
log(z) = RHOZ*log(z(-1)) + SIGZ/100*epsz;

@#if PREFSHOCK == 1
[name='Preference shifter process']
log(zeta) = RHOZETA*log(zeta(-1)) + SIGZETA/100*epszeta;
@#endif

[name='Output growth (q-on-q)']
YGR = GAMQ + 100*(log(y/steady_state(y)) - log(y(-1)/steady_state(y)) + log(z/steady_state(z)))
@#if MEASERR == 1
    + SIGYGR*measerr_YGR
@#endif
;

[name='Annualized inflation']
INFL = PA + 400*log(p/steady_state(p))
@#if MEASERR == 1
    + SIGINFL*measerr_INFL
@#endif
;

[name='Annualized nominal interest rate']
INT = PA + RA + 4*GAMQ + 400*log(r/steady_state(r))
@#if MEASERR == 1
    + SIGINT*measerr_INT
@#endif
;

end; % [model] end


% =========================================================================
% Steady state Model
% =========================================================================
steady_state_model;
GAMMA    = 1+GAMQ/100;
BETA     = 1/(1+RA/400);
PBARSTAR = 1+PA/400;
@#if PREFSHOCK==1
ZETABAR = 1;
zeta = ZETABAR;
@#endif
ZBAR = 1;
z    = ZBAR;
p    = PBARSTAR;
g    = 1/C_o_Y;
r    = GAMMA*ZBAR*PBARSTAR/BETA;
c    = (1-NU)^(1/TAU);
y    = g*c;
YGR  = GAMQ;
INFL = PA;
INT  = PA + RA + 4*GAMQ;
end; % [steady_state_model] end


% =========================================================================
% Declare settings for shocks
% =========================================================================
shocks;
var epsr = 1;
var epsg = 1;
var epsz = 1;
@#if PREFSHOCK == 1
var epszeta = 1;
@#endif
@#if MEASERR == 1
var measerr_YGR = 1;
var measerr_INFL = 1;
var measerr_INT = 1;
@#endif
end; % [shocks] end


% =========================================================================
% Specify Priors
% =========================================================================
estimated_params;
% --------------------------------------------------------------------------------------------------------------------------------------------------
%PARAMETER_NAME, INITIAL_VALUE, LOWER_BOUND, UPPER_BOUND, PRIOR_SHAPE,   PRIOR_MEAN, PRIOR_STANDARD_ERROR, PRIOR_3RD_PARAMETER, PRIOR_4TH_PARAMETER;
% --------------------------------------------------------------------------------------------------------------------------------------------------
RA,              1,             1e-5,        10,          gamma_pdf,     0.8,        0.5;
PA,              3.2,           1e-5,        20,          gamma_pdf,     4,          2;
GAMQ,            0.55,          -5,          5,           normal_pdf,    0.4,        0.2;
TAU,             2,             1e-5,        10,          gamma_pdf,     2,          0.5;
NU,              0.1,           1e-5,        0.99999,     beta_pdf,      0.1,        0.05;
PSIP,            1.5,           1e-5,        10,          gamma_pdf,     1.5,        0.25;
PSIY,            0.125,         1e-5,        10,          gamma_pdf,     0.5,        0.25;
@#if MONPOL == 3
PSIDY,           0.2,           1e-5,        10,          gamma_pdf,     0.2,        0.15;
@#endif
RHOR,            0.75,          1e-5,        0.99999,     beta_pdf,      0.5,        0.2;
RHOG,            0.95,          1e-5,        0.99999,     beta_pdf,      0.8,        0.1;
RHOZ,            0.9,           1e-5,        0.99999,     beta_pdf,      0.66,       0.15;
SIGR,            0.2,           1e-8,        5,           inv_gamma_pdf, 0.3,        4;
SIGG,            0.6,           1e-8,        5,           inv_gamma_pdf, 0.4,        4;
SIGZ,            0.3,           1e-8,        5,           inv_gamma_pdf, 0.4,        4;
@#if INDEXATION == 1
IOTAP,           0.5,           1e-8,        0.99999,     beta_pdf,      0.5,        0.15;
@#endif
@#if PREFSHOCK == 1
RHOZETA,         0.75,          1e-5,        0.99999,     beta_pdf,      0.5,        0.2;
SIGZETA,         0.2,           1e-8,        5,           inv_gamma_pdf, 0.3,        4;
@#endif
@#if MEASERR == 1
SIGYGR,          0.23,          1e-8,        5,           inv_gamma_pdf, 0.3,        4;
SIGINFL,         0.56,          1e-8,        5,           inv_gamma_pdf, 0.3,        4;
SIGINT,          0.66,          1e-8,        5,           inv_gamma_pdf, 0.3,        4;
@#endif
end; % [estimated_params] end


@#if COMPUTATIONS >= 0
% =========================================================================
% Steady-State, Checks and Diagnostics
% =========================================================================
steady;                 % compute steady state given the starting values
resid;                  % check the residuals of model equations evaluated at steady state
check;                  % check Blanchard-Kahn conditions
model_diagnostics;      % check obvious model errors
disp(' ');
@#endif

@#if COMPUTATIONS == 1
% =========================================================================
% Lack of Identification Analysis for all possible combinations of VAROBS
% =========================================================================
RESULTS = {};
TABLE = {};
irun = 0;
for iset = 1:@{VAROBSCOMBINATIONS}
    varobsset = nchoosek(1:M_.endo_nbr,iset);
    varobsset_nbr = size(varobsset,1);
    RESULTS{iset}=cell(varobsset_nbr,1);
    for ii = 1:varobsset_nbr
        irun = irun+1;
        options_.varobs = M_.endo_names(varobsset(ii,:),:); % set VAROBS
        fprintf('**** Null Space VAROBS: %s ****' , strjoin(options_.varobs));
        identification(ar=30,no_identification_strength, parameter_set=calibration); % run identification analysis
        fprintf('**** Bruteforce VAROBS: %s ****' , strjoin(options_.varobs));
        identification(ar=30,no_identification_strength, parameter_set=calibration,checks_via_subsets=1, max_dim_subsets_groups=20); % run identification analysis
        fprintf('**** Monte Carlo Testing VAROBS: %s ****' , strjoin(options_.varobs));
        identification(ar=30,no_identification_strength, prior_mc=100, nograph, no_identification_minimal); % run identification analysis        
        % save results into structure
        RESULTS{iset}{ii} = load([M_.fname, '/identification/', M_.fname, '_calibration_identif.mat'], 'ide_moments_point', 'ide_spectrum_point', 'ide_minimal_point');
        RESULTS{iset}{ii}.varobs = options_.varobs;
        
        % find out if PSIP, PSIY, RHOR, SIGR are among problematic parameter sets
        idx_PSIP = find(estim_params_.param_vals(:,1)==find(contains(M_.param_names,'PSIP')));
        idx_PSIY = find(estim_params_.param_vals(:,1)==find(contains(M_.param_names,'PSIY')));
        idx_RHOR = find(estim_params_.param_vals(:,1)==find(contains(M_.param_names,'RHOR')));
        idx_SIGR = find(estim_params_.param_vals(:,1)==find(contains(M_.param_names,'SIGR')));
        for jide = 1:3 %loop over moments, minimal and spectrum identification results
            if jide == 1
                ide = RESULTS{iset}{ii}.ide_moments_point;
            elseif jide == 2
                ide = RESULTS{iset}{ii}.ide_minimal_point;
            elseif jide == 3
                ide = RESULTS{iset}{ii}.ide_spectrum_point;
            end
            
            if isfield(ide,'rank') %check if anything went wrong
                indx_sets = find(~cellfun(@isempty,ide.problpars));
                if isempty(indx_sets)
                    TABLE{irun,jide} = '$\checkmark\checkmark$'; %all parameters are identified
                else
                    if ~any(cellfun(@(x) nnz(ismember(x,[idx_PSIP, idx_PSIY, idx_RHOR, idx_SIGR])),ide.problpars(indx_sets), 'UniformOutput', 1))
                        TABLE{irun,jide} = '$\checkmark$'; % PSIP,PSIY,RHOR,SIGR are identified, but other parameters are not
                    else
                        TABLE{irun,jide} = '$[';
                        if any(cellfun(@(x) nnz(ismember(x,idx_PSIP)),ide.problpars(indx_sets), 'UniformOutput', 1))
                            TABLE{irun,jide} = [TABLE{irun,jide} '\psi_\pi ']; % PSIP is within unidentified sets
                        end
                        if any(cellfun(@(x) nnz(ismember(x,idx_PSIY)),ide.problpars(indx_sets), 'UniformOutput', 1))
                            TABLE{irun,jide} = [TABLE{irun,jide} '\psi_y ']; % PSIY is within unidentified sets
                        end
                        if any(cellfun(@(x) nnz(ismember(x,idx_RHOR)),ide.problpars(indx_sets), 'UniformOutput', 1))
                            TABLE{irun,jide} = [TABLE{irun,jide} '\rho_R ']; % RHOR is within unidentified sets
                        end
                        if any(cellfun(@(x) nnz(ismember(x,idx_SIGR)),ide.problpars(indx_sets), 'UniformOutput', 1))
                            TABLE{irun,jide} = [TABLE{irun,jide} '\sigma_R ']; % SIGR is within unidentified sets
                        end
                        TABLE{irun,jide} = [TABLE{irun,jide} ']$'];
                    end
                end
            else
                TABLE{irun,jide} = 'err';
            end
            TABLE{irun,4} = ['$', strjoin(M_.endo_names_tex(varobsset(ii,:)),','), '$'];
        end
    end
end
TABLE = cell2table(TABLE,'VariableNames',{'Moments','Minimal','Spectrum','Varobs'});
save('RESULTS_TABLE.mat','RESULTS', 'TABLE');  % save results as table

% Now use this matlab table as input for (slightly modified) latexTable.m
[~,pname] = fileparts(pwd); %folder name
inputTbl.data = TABLE;
inputTbl.dataFormat = {'%s'};
inputTbl.tableColumnAlignment = 'c';
inputTbl.tableCaption = strrep(pname,'_',' ');
inputTbl.tableBorders = 1;    
inputTbl.longtable = 1;
inputTbl.makeCompleteLatexDocument = 1;    
latex = LatexTable(inputTbl);
filenam = sprintf('table.tex'); 

% save LaTex code as file and compile it
fid=fopen(filenam,'w');
[nrows,ncols] = size(latex);
for row = 1:nrows
    fprintf(fid,'%s\n',latex{row,:});
end
fclose(fid);
system(['pdflatex -halt-on-error ' filenam])

@#endif


@#if COMPUTATIONS == 2
% =========================================================================
% Simulate data
% =========================================================================
stoch_simul(order=1,IRF=0,periods=@{SAMPLESIZE});
save('../simdat.mat','YGR','INFL','INT','y','c','r','p','g','z');
@#if PREFSHOCK == 1
    save('../simdat.mat','zeta','-append');
@#endif
@#endif


@#if COMPUTATIONS == 3
% =========================================================================
% Bayesian Learning Rate Indicator
% =========================================================================
varobs @{VAROBS} ;
copyfile('../simdat.mat');
pause(1);

@#if HESSIAN == 1 && FINDMODEADVANCED == 1
optimtypevals=[9 8 4 7 1 6]; % Declare types of Dynare optimisers used

for jj=1:length(optimtypevals)
    numopt = optimtypevals(jj);
    if     optimtypevals(jj)==1, disp('RUNNING OPTIMISER-TYPE: Matlab fmincon');
    elseif optimtypevals(jj)==2, disp('RUNNING OPTIMISER-TYPE: Simulated annealing');
    elseif optimtypevals(jj)==4, disp('RUNNING OPTIMISER-TYPE: Sims');
    elseif optimtypevals(jj)==5, disp('RUNNING OPTIMISER-TYPE: Marco Ratto optimiser');
    elseif optimtypevals(jj)==6, disp('RUNNING OPTIMISER-TYPE: Dynare Monte Carlo');
    elseif optimtypevals(jj)==7, disp('RUNNING OPTIMISER-TYPE: fminsearch');
    elseif optimtypevals(jj)==8, disp('RUNNING OPTIMISER-TYPE: Nelson-Mead Simplex');
    elseif optimtypevals(jj)==9, disp('RUNNING OPTIMISER-TYPE: CMA-ES');
    end

    if jj==1
        eval(['options_.mode_compute=' num2str(numopt) ';']) ;
        estimation(cova_compute=1,mh_replic=0,plot_priors=0,datafile=simdat,first_obs=@{FIRSTOBS},nobs=@{SAMPLESIZE},tex);
        movefile([M_.fname '_mode.mat'], 'newmode.mat'); % rename mode file
    elseif jj>1 && jj<length(optimtypevals)
        eval(['options_.mode_compute=' num2str(numopt) ';']) ;
        estimation(mode_file=newmode,cova_compute=1,mh_replic=0,plot_priors=0,datafile=simdat,first_obs=@{FIRSTOBS},nobs=@{SAMPLESIZE},tex);
        movefile([M_.fname '_mode.mat'], 'newmode.mat'); % rename mode file
    elseif jj==length(optimtypevals) && numopt == 6
        eval(['options_.mode_compute=' num2str(numopt) ';']) ;
        options_.gmhmaxlik.target=0.2;
        options_.gmhmaxlik.iterations = 1;
        estimation(optim=('AcceptanceRateTarget',0.2),mode_check,mode_file=newmode,cova_compute=1,mh_replic=0,plot_priors=0,datafile=simdat,first_obs=@{FIRSTOBS},nobs=@{SAMPLESIZE},tex);
    end
end
movefile([M_.fname '_mode.mat'], 'finalmode.mat'); % rename final mode file
% Save Bayesian Learning Rate Indicator based on the Hessian. Note that if FINDMODEADVANCED=1, 
% the Hessian is found using a Monte-Carlo based optimization routine (mode compute 6)
% This is in general not a good approximation of the true Hessian, but is
% useful for the MCMC estimation. See https://git.dynare.org/Dynare/dynare/wikis/mode-compute-6
weakresults_hessian = (diag(oo_.posterior.optimization.Variance).*@{SAMPLESIZE}).^(-1)'; 
@#endif


@#if HESSIAN == 1 && FINDMODEADVANCED == 0
estimation(datafile=simdat,
           first_obs=@{FIRSTOBS},
           nobs=@{SAMPLESIZE},
           mode_compute=@{OPTIMIZER},
           mode_check,
           plot_priors=0,
           mh_replic=0,
           cova_compute=1,
           tex);
movefile([M_.fname '_mode.mat'], 'finalmode.mat'); % rename mode file
% Save Bayesian Learning Rate Indicator based on the Hessian
weakresults_hessian = (diag(oo_.posterior.optimization.Variance).*@{SAMPLESIZE}).^(-1)';
@#endif


@#if MCMC == 1
estimation(datafile=simdat,
           first_obs=@{FIRSTOBS},
           nobs=@{SAMPLESIZE},
           mode_compute=0, % do not compute mode again
           mode_file=finalmode, % use already computed mode
           plot_priors=0,
           mh_replic=@{MH_REPLIC},
           mh_nblocks=@{MH_NBLOCKS},
           mh_jscale=@{JSCALEINTEGER}.@{JSCALEDECIMAL}, % jscale is set to JSCALEINTEGER.JSCALEDECIMAL such that average acceptance ratios lie between 0.2 and 0.3
           tex);
% Save Bayesian Learning Rate Indicator based on MCMC
weakresults_mcmc = (diag(oo_.posterior.metropolis.Variance).*@{SAMPLESIZE}).^(-1)';
@#endif

collect_latex_files;
system(['pdflatex -halt-on-error -interaction=batchmode ' M_.fname '_TeX_binder.tex'])

@#endif
