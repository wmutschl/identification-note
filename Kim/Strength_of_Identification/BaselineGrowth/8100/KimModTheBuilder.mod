% Replication files for:
% Ivashchenko and Mutschler - On the effect of observables, functional specifications, 
% and modeling choice on local identification in linearized DSGE models
% =========================================================================
% Note that this Mod file needs to be run with at least Dynare 4.6-unstable
% =========================================================================
% This mod file implements 
% (1) several different features inside the mode of Kim (2003) - Functional
% equivalence between intertemporal and multisectoral investment adjustment
% costs; namely:
%   - Investment adjustment costs: none, multisectoral, intertemporal (growth or level), or both types
%   - utility function: Logarithmic or CRRA
%   - consumption habit: none, external or internal
%   - utiliation costs of capital: with or without
%   - investment specific technological shock: with or without
%   - labor: with or without
%   - monetary policy: with or without
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
% - MSADJ: Multisectoral adjustment cost specification
%   - 0: no multisectoral adjustment costs, i.e. Y_t^d = C_t + I_t + \Psi^K_t(U^K_t) K_{t-1}
%   - 1: with multisectoral adjustment costs, i.e. Y_t^d = ((1-SAV)*(C_t/(1-SAV))^(1-\theta) + SAV*(I_t/SAV)^(1+\theta) )^(1/(1+\theta)) + \Psi^K_t(U^K_t) K_{t-1} with SAV = I/Y
% -------------------------------------------------------------------------
% - IAC: Investment adjustment cost specification
%   - 0: no investment adjustment costs, i.e. K_t = (1-\delta) K_{t-1} + \upsilon_t I_t
%   - 1: investment level costs a la KIM (2003), i.e. K_t = ( (1-\delta) K_{t-1}^(1-\kappa) + \delta (\upsilon_t I_t/\delta)^(1-\kappa) )^(1/(1-\kappa))
%   - 2: investment growth costs a la CEE (2005), i.e. K_t = (1-\delta) K_{t-1} + \upsilon_t I_t (1-S_t), where S_t is a function of investment growth (I_t/I_{t-1})
% -------------------------------------------------------------------------
% - UTIL: Functional form of utility function
%   - 0: LOGLOG, i.e. U(C_t-H_t)=\log(C_t-H_t)+\gamma_L \log(1-L_t)
%   - 1: CRRA, i.e. U(C_t-H_t)=(C_t-H_t)^(1-\eta_C)/(1-\eta_C)+\gamma_L (1-L_t)^(1-\eta_L)/(1-\eta_L)
% -------------------------------------------------------------------------
% - HABIT: Functional form of habit formation
%   - 0: no internal habit, i.e. H_t = 0
%   - 1: with external habit, i.e. H_t = \hbar C
%   - 2: with internal habit, i.e. H_t = \hbar C_{t-1}
% -------------------------------------------------------------------------
% - CAPUTIL: Model with or without utiliation costs of capital
%   - 0: no capital utilization, i.e. \Psi^K_t(U_t^K) = 0
%   - 1: with capital utilization, i.e. \Psi^K_t(U_t^K) = (1-\psi^K)(U_t^K-U^K) + \psi^K/2(U^K_t-U^K)^2
% -------------------------------------------------------------------------
% - INVESTSHOCK: Model with or without investment specific technological shock
%   - 0: no investment specific technological shock, i.e. \upsilon_t = 1
%   - 1: with investment specific technological shock, i.e. \log{\upsilon_t} = \rho_\upsilon \log{\upsilon_{t-1}} + \sigma_\upsilon \varepsilon_t^\upsilon
% -------------------------------------------------------------------------
% - LABOR: Model with or without labor choice
%   - 0: no labor choice
%   - 1: with labor choice
% -------------------------------------------------------------------------
% - MONPOL: Specification of Taylor Rule
%   - 0: no intertemporal bond choice and no monetary policy
%   - 1: intertemporal bond choice and monetary policy that reacts to deviations in output growth and inflation target
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
y        ${Y}$        (long_name='output')
c        ${C}$        (long_name='consumption')
iv       ${I}$        (long_name='investment')
rk       ${R^{K}}$    (long_name='rental rate of capital')
k        ${K}$        (long_name='private capital stock')
lam      ${\Lambda}$  (long_name='marginal utility, i.e. Lagrange multiplier budget')
q        ${Q}$        (long_name='Tobins Q, i.e. Lagrange multiplier capital stock')
a        ${A}$        (long_name='total factor productivity')
@#if CAPUTIL == 1
uk       ${U^K}$      (long_name='capital utilization rate')
@#endif
@#if INVESTSHOCK == 1
upsilon  ${\upsilon}$ (long_name='investement-specific technology')
@#endif
@#if LABOR == 1
l        ${L}$        (long_name='labor')
w        ${W}$        (long_name='real wage')
@#endif
@#if MONPOL == 1
r        ${R}$        (long_name='nominal interest rate')
pie      ${\pi}$     (long_name='inflation rate')
@#endif
; % [var] end


% =========================================================================
% Declare exogenous variables
% =========================================================================
varexo
epsa       ${\varepsilon^A}$        (long_name='total factor productivity shock')
@#if INVESTSHOCK == 1
epsupsilon ${\varepsilon^\upsilon}$ (long_name='investement-specific technology shock')
@#endif
@#if MONPOL == 1
epsr       ${\varepsilon^R}$        (long_name='monetary policy shock')
@#endif
; % [varexo] end


% =========================================================================
% Declare parameters
% =========================================================================
parameters
ALPHA      ${\alpha}$           (long_name='bias towards capital in production')
RA         ${r_{A}}$            (long_name='annual steady-state real interest rate (defines discount factor)')
DELTA      ${\delta}$           (long_name='depreciation rate')
QBAR       ${\bar{Q}}$          (long_name='steady state Tobins Q')
ABAR       ${\bar{A}}$          (long_name='steady state technology')
RHOA       ${\rho_A}$           (long_name='persistence TFP')
SIGA       ${\sigma_A}$         (long_name='standard deviation TFP shock')
@#if MSADJ == 1
THETA      ${\theta}$           (long_name='multisectoral adjustment cost parameter')
@#endif
@#if IAC > 0
KAPPA      ${\kappa}$           (long_name='investment adjustment cost parameter')
@#endif
@#if HABIT > 0
HBAR       ${\hbar}$            (long_name='consumption habit parameter')
@#endif
@#if UTIL == 1
ETAC       ${\eta^c}$           (long_name='consumption utility parameter')
@#endif
@#if CAPUTIL == 1
UKBAR      ${\bar{u}^K}$        (long_name='steady state utilization rate')
PSIK       ${\psi^K}$           (long_name='capital utilization parameter')
@#endif
@#if INVESTSHOCK == 1
UPSILONBAR ${\bar{\upsilon}}$   (long_name='steady state investment specific technology')
RHOUPSILON ${\rho_\upsilon}$    (long_name='persistence investment specific technology shock')
SIGUPSILON ${\sigma_\upsilon}$  (long_name='standard deviation investment specific technology shock')
@#endif
@#if LABOR == 1
GAMMAL     ${\gamma_L}$         (long_name='labor weight in utility')
@#endif
@#if LABOR == 1 && UTIL == 1
ETAL       ${\eta_L}$           (long_name='labor disutility parameter')
@#endif
@#if MONPOL == 1
PIEA       ${\pi_{A}}$          (long_name='target inflation rate')
PSIPIE     ${\psi_\pi}$         (long_name='Taylor rule inflation sensitivity')
PSIY       ${\psi_Y}$           (long_name='Taylor rule output sensitivity')
RHOR       ${\rho_R}$           (long_name='persistence Taylor rule')
SIGR       ${\sigma_R}$         (long_name='standard deviation monetary policy shock')
@#endif
; % [parameters] end


% =========================================================================
% Calibrate parameter values and compute implied steady state (denoted here with BARS)
% =========================================================================
% Relevant calibration for full model
IBAR_O_YBAR = 0.25;         % average investment output ratio I/Y
KBAR_O_YBAR = 10;           % average capital capital productivity K/Y
DELTA       = IBAR_O_YBAR/KBAR_O_YBAR;  % quarterly depreciation rate
RA          = 2;            % annual nominal interest rate
BETTABAR    = 1/(1+RA/400); % discount factor 
THETA       = 1.5;          % multisectoral investment adjustment cost parameter
KAPPA       = 2;            % intertemporal investment adjustment cost parameter
HBAR        = 0.6;          % consumption habit parameter
ETAC        = 2;            % risk aversion parameter of households
ETAL        = 1.5;          % inverse of intertemporal elasticity of substitution w.r.t leisure
@#if LABOR == 1
LBAR        = 0.33;         % steady state labor
@#else
LBAR        = 1;            % normalize labor
@#endif
RHOA        = 0.5;          % technology persistence
RHOUPSILON  = 0.5;          % investment-specific technological progress persistence
RHOR        = 0.5;          % persistence of Taylor rule
SIGA        = 0.6;          % technology standard deviation
SIGUPSILON  = 0.6;          % investment-specific technological progress standard deviation
SIGR        = 0.2;          % standard deviation of monetary policy shock
PIEA        = 3.2;          % annual inflation target
PIEBAR      = 1+PIEA/400;   % central bank target
PSIPIE      = 1.5;          % inflation sensitivity of Taylor rule
PSIY        = 0.125;        % output gap sensitivity of Taylor rule
ABAR        = 1;            % normalize steady state technology level
QBAR        = 1;            % normalize steady state Tobin's Q
UKBAR       = 1;            % normalize steady state capital utilization rate
UPSILONBAR  = 1;            % normalize investment-specific technology
@#if CAPUTIL == 1
PSIK = ((UKBAR^2-3)/(2*UKBAR))^(-1)*( (1/BETTABAR+DELTA-1)*QBAR/UKBAR-1/UKBAR); % PSIK is implicitly defined via steady state from foc wrt UK and FOC wrt K
PSSIKBAR = (1-PSIK)*(UKBAR-UKBAR) + PSIK/2*(UKBAR-UKBAR)^2;                     % this is just 0
RKBAR = (1/BETTABAR+DELTA-1)*QBAR/UKBAR + PSSIKBAR/UKBAR;                       % foc K
@#else
RKBAR = (1/BETTABAR+DELTA-1)*QBAR; % foc K
@#endif
ALPHA = RKBAR*KBAR_O_YBAR;                         % capital share in production
KSBAR_o_LBAR = (ALPHA*ABAR/RKBAR)^(1/(1-ALPHA));   % steady state Ks/L from capital demand
KSBAR = KSBAR_o_LBAR*LBAR;                         % steady state Ks from identity
@#if CAPUTIL == 1
KBAR_o_LBAR = KSBAR_o_LBAR/UKBAR; % utilized K/L
KBAR = KSBAR/UKBAR;               % utilized K
@#else
KBAR_o_LBAR = KSBAR_o_LBAR;       % utilized K/L
KBAR = KSBAR;                     % utilized K
@#endif
IVBAR = DELTA/UPSILONBAR*KBAR;          % steady state I from capital accumulation
YBAR = ABAR*KSBAR^ALPHA*LBAR^(1-ALPHA); % steady state Y from production function
@#if CAPUTIL == 1
YSBAR = YBAR-PSSIKBAR*KBAR;     % steady state C from market clearing
@#else
YSBAR = YBAR;                   % steady state C from market clearing
@#endif
SAVBAR = IVBAR / YSBAR;  % definition of SAV
CBAR = (1-SAVBAR)*YSBAR; % identity from market clearing
@#if MSADJ == 1
YDBAR = ( (1-SAVBAR)*(CBAR/(1-SAVBAR))^(1+THETA) + SAVBAR*(IVBAR/SAVBAR)^(1+THETA) )^(1/(1+THETA)); % demand side of budget restriction
XCBAR = ( CBAR/((1-SAVBAR)*YDBAR) )^THETA;                                                          % auxiliary expression
XIVBAR = ( IVBAR/(SAVBAR*YDBAR) )^THETA;                                                            % auxiliary expression
@#else
YDBAR = CBAR + IVBAR; % demand side of budget restriction
XCBAR = 1;            % auxiliary expression
XIVBAR = 1;           % auxiliary expression
@#endif
@#if HABIT == 0
HABBAR = 0;          % no habit
DHABP_DC_BAR = 0;    % auxiliary expression
@#endif
@#if HABIT == 1
HABBAR = HBAR*CBAR;  % external habit in steady state
DHABP_DC_BAR = 0;    % auxiliary expression
@#endif
@#if HABIT == 2
HABBAR = HBAR*CBAR;  % internal habit in steady state
DHABP_DC_BAR = HBAR; % auxiliary expression
@#endif
@#if UTIL == 0
DU_DC_BAR = (CBAR-HABBAR)^(-1);  % auxiliary derivative \partial U_t / \partial C_t in steady state
DUP_DC_BAR = (CBAR-HABBAR)^(-1)*(-DHABP_DC_BAR); % auxiliary derivative \partial U_{t+1} / \partial C_t in steady state
@#endif
@#if UTIL == 1
DU_DC_BAR = (CBAR-HABBAR)^(-ETAC);  % auxiliary derivative \partial U_t / \partial C_t in steady state
DUP_DC_BAR = (CBAR-HABBAR)^(-ETAC)*(-DHABP_DC_BAR); % auxiliary derivative \partial U_{t+1} / \partial C_t in steady state
@#endif
LAMBAR = XCBAR^(-1)*(DU_DC_BAR + BETTABAR*DUP_DC_BAR); % foc C
@#if LABOR == 1
WBAR = (1-ALPHA)*ABAR*(KSBAR_o_LBAR)^ALPHA; % labor demand
    @#if UTIL == 0
GAMMAL = LAMBAR*WBAR*(1-LBAR);      % foc L
    @#endif
    @#if UTIL == 1
GAMMAL = LAMBAR*WBAR*(1-LBAR)^ETAL; % foc L
    @#endif
@#endif
@#if MONPOL == 1
RBAR = PIEBAR/BETTABAR; % Steady state interest rate
@#endif


% =========================================================================
% Model equations
% =========================================================================
model;
% -------------------------------------------------------------------------
% Auxiliary parameters
#BETA = 1/(1+RA/400);
% -------------------------------------------------------------------------
% Auxiliary expression for investment specific technology shock
@#if INVESTSHOCK == 0
#upsilon = 1;
#upsilonp = 1;
@#else
#upsilonp = upsilon(+1);
@#endif
% -------------------------------------------------------------------------
% Auxiliary expressions for utilized capital
@#if CAPUTIL == 1
#ks = uk(+1)*k;
#ksback = uk*k(-1);
#Psik = (1-PSIK)*(uk-steady_state(uk)) + PSIK/2*(uk-steady_state(uk))^2;
#Psikp = (1-PSIK)*(uk(+1)-steady_state(uk)) + PSIK/2*(uk(+1)-steady_state(uk))^2;
#Psikprime = (1-PSIK) + PSIK*(uk-steady_state(uk));
#Psikprimep = (1-PSIK) + PSIK*(uk(+1)-steady_state(uk));
#Psikbar = (1-PSIK)*(steady_state(uk)-steady_state(uk)) + PSIK/2*(steady_state(uk)-steady_state(uk))^2;
@#else
#ks = k;
#ksback = k(-1);
@#endif
% -------------------------------------------------------------------------
% Auxiliary expressions for intertemporal investment adjustment cost
@#if IAC == 2
#S = KAPPA/2*(iv/iv(-1)-1)^2;
#Sprime = KAPPA*(iv/iv(-1)-1);
#Sprimep = KAPPA*(iv(+1)/iv-1);
#Sprimeprime = KAPPA;
@#endif
% -------------------------------------------------------------------------
% Auxiliary expressions for multisectoral adjustment cost
@#if MSADJ == 1
#SAVBAR = steady_state(iv)/( steady_state(y)
	@#if CAPUTIL == 1
    -Psikbar*steady_state(k)
    @#endif
);
#yd = ((1-SAVBAR)*(c/(1-SAVBAR))^(1+THETA) + SAVBAR*(iv/SAVBAR)^(1+THETA))^(1/(1+THETA));
#xc = ( c/((1-SAVBAR)*yd) )^THETA;
#xiv = ( iv/(SAVBAR*yd) )^THETA;
@#else
#yd = c + iv;
#xc = 1;
#xiv = 1;
@#endif
% -------------------------------------------------------------------------
% Auxiliary expressions for no habit
@#if HABIT == 0
#hab = 0;
#habp = 0;
#dhab_dc = 0;
#dhabp_dc = 0;
@#endif
% Auxiliary expressions for external habit
@#if HABIT == 1
#hab = HBAR*steady_state(c);
#habp = HBAR*steady_state(c);
#dhab_dc = 0;
#dhabp_dc = 0;
@#endif
% Auxiliary expressions for internal habit
@#if HABIT == 2
#hab = HBAR*c(-1);
#habp = HBAR*c;
#dhab_dc = 0;
#dhabp_dc = HBAR;
@#endif
% -------------------------------------------------------------------------
% Auxiliary expressions for log utility function
@#if UTIL == 0
#du_dc = (c-hab)^(-1)*(1-dhab_dc);
#dup_dc = (c(+1)-habp)^(-1)*(-dhabp_dc);
    @#if LABOR == 1
#du_dl = GAMMAL*(1-l)^(-1)*(-1);
    @#endif
@#endif
% Auxiliary expressions for CRRA utility function
@#if UTIL == 1
#du_dc = (c-hab)^(-ETAC)*(1-dhab_dc);
#dup_dc = (c(+1)-habp)^(-ETAC)*(-dhabp_dc);
    @#if LABOR == 1
#du_dl = GAMMAL*(1-l)^(-ETAL)*(-1);
    @#endif
@#endif
% -------------------------------------------------------------------------
% Auxiliary expressions for production function
@#if LABOR == 0
#pf = ksback^ALPHA;
#dpf_dksback = ALPHA*ksback^(ALPHA-1);
@#endif
@#if LABOR == 1
#pf = ksback^ALPHA*l^(1-ALPHA);
#dpf_dl = (1-ALPHA)*(ksback/l)^ALPHA;
#dpf_dksback = ALPHA*(ksback/l)^(ALPHA-1);
@#endif
% -------------------------------------------------------------------------
% Auxiliary expressions for monetary policy rule
@#if MONPOL == 1
#Piebar = 1+PIEA/400;
#rstar = steady_state(r)*(pie/Piebar)^PSIPIE*(y/steady_state(y))^PSIY;
@#endif
% -------------------------------------------------------------------------


[name='foc household wrt c (marginal utility of consumption)']
xc*lam = du_dc + BETA*dup_dc;

[name='foc household wrt iv (Tobins Q)']
xiv*lam = lam*q*upsilon
    @#if IAC == 1
    ^(-KAPPA)*(DELTA*k/iv)^KAPPA
    @#endif
    @#if IAC == 2
    *(1-S-(iv/iv(-1))*Sprime) + BETA*lam(+1)*q(+1)*upsilonp*(iv(+1)/iv)^2*Sprimep
    @#endif
;

[name='foc household wrt k (Euler equation capital)']
lam*q = BETA*lam(+1)*(rk(+1)
    @#if CAPUTIL == 1
	*uk(+1) - Psikp
    @#endif
        + (1-DELTA)*q(+1)
            @#if IAC == 1
            *(k(+1)/k)^KAPPA
            @#endif
);

@#if CAPUTIL == 1
[name='foc household wrt uk (optimal capital utilization)']
rk = Psikprime;
@#endif

[name='capital accumulation']
@#if IAC == 0
k = (1-DELTA)*k(-1) + upsilon*iv;
@#endif
@#if IAC == 1
k = ((1-DELTA)*k(-1)^(1-KAPPA) + DELTA*(upsilon*iv/DELTA)^(1-KAPPA))^(1/(1-KAPPA));
@#endif
@#if IAC == 2
k = (1-DELTA)*k(-1) + upsilon*iv*(1-S);
@#endif

[name='production function']
y = a*pf;

[name='firm capital demand']
rk = a*dpf_dksback;

[name='market clearing (resource constraint)']
y = yd
    @#if CAPUTIL == 1
	+ Psik*k(-1)
    @#endif
;

[name='Evolution of technology']
log(a) = RHOA*log(a(-1)) + SIGA*epsa;

@#if LABOR == 1
[name='foc household wrt l']
lam*w = -du_dl;
@#endif

@#if LABOR == 1
[name='labor demand']
w = a*dpf_dl;
@#endif

@#if INVESTSHOCK == 1
[name='Evolution of investment specific technology shock']
log(upsilon) = RHOUPSILON*log(upsilon(-1)) + SIGUPSILON*epsupsilon;
@#endif

@#if MONPOL == 1
[name='foc household wrt b (Euler equation bonds)']
lam = BETA*r/pie(+1)*lam(+1);
@#endif

@#if MONPOL == 1
[name='Taylor Rule']
r = rstar^(1-RHOR)*r(-1)^RHOR*exp(SIGR*epsr);
@#endif

end; % [model] end


% =========================================================================
% Steady State Model
% =========================================================================
@#if UTIL == 0 || (UTIL==1 && LABOR == 0)
steady_state_model;
BETTA = 1/(1+RA/400);
A = ABAR;
a = A;
Q = QBAR;
q = Q;
% -------------------------------------------------------------------------
@#if INVESTSHOCK == 1
UPSILON = UPSILONBAR;
upsilon = 1;
@#else
UPSILON = 1;
@#endif
% -------------------------------------------------------------------------
@#if CAPUTIL == 1
RK = 1-PSIK;
UK = (1/BETTA+DELTA-1)*Q/RK;
uk = UK;
rk = RK;
PSSIK = (1-PSIK)*(UK-UK) + PSIK/2*(UK-UK)^2;
@#else
RK = (1/BETTA+DELTA-1)*Q; % foc k
@#endif
rk = RK;
% -------------------------------------------------------------------------
KS_o_L = ((ALPHA*A)/RK)^(1/(1-ALPHA)); % capital demand
@#if CAPUTIL == 1
K_o_L = KS_o_L/UK; % Utilized capital
@#else
K_o_L = KS_o_L; % Utilized capital
@#endif
% -------------------------------------------------------------------------
Y_o_L = A*KS_o_L^ALPHA;       % production function
@#if CAPUTIL == 1
YS_o_L = Y_o_L - PSSIK*K_o_L; % market clearing
@#else
YS_o_L = Y_o_L; % market clearing
@#endif
% -------------------------------------------------------------------------
IV_o_L = DELTA*K_o_L/UPSILON; % capital accumulation
SAV = IV_o_L / YS_o_L;
C_o_L = (1-SAV)*YS_o_L;
% -------------------------------------------------------------------------
@#if MSADJ == 1
YD_o_L = ( (1-SAV)*(C_o_L/(1-SAV))^(1+THETA) + SAV*(IV_o_L/SAV)^(1+THETA) )^(1/(1+THETA));
XC = ( C_o_L/((1-SAV)*YD_o_L) )^THETA;
XIV = ( IV_o_L/(SAV*YD_o_L) )^THETA;
@#else
YD_o_L = C_o_L + IV_o_L;
XC = 1;
XIV = 1;
@#endif
% -------------------------------------------------------------------------
@#if HABIT == 0
HAB_o_L = 0;           % no habit
DHABP_DC = 0;          % auxiliary expression
@#endif
@#if HABIT == 1
HAB_o_L = HBAR*C_o_L;  % external habit in steady state
DHABP_DC = 0;          % auxiliary expression
@#endif
@#if HABIT == 2
HAB_o_L = HBAR*C_o_L;  % internal habit in steady state
DHABP_DC = HBAR;       % auxiliary expression
@#endif
% -------------------------------------------------------------------------
% The labor level in closed-form for UTIL=0
@#if LABOR == 0
L = 1;
@#else
W = (1-ALPHA)*A*KS_o_L^ALPHA; % labor demand
w = W;
@#if UTIL == 0
L = XC^(-1)/GAMMAL*(C_o_L-HAB_o_L)^(-1)*W*(1-BETTA*DHABP_DC) /( 1 + XC^(-1)/GAMMAL*(C_o_L-HAB_o_L)^(-1)*W*(1-BETTA*DHABP_DC) );
@#endif
l = L;
@#endif
%--------------------------------------------------------------------------
@#if UTIL == 0
DU_DC = (C_o_L*L-HAB_o_L*L)^(-1);                 % auxiliary derivative \partial U_t / \partial C_t in steady state
DUP_DC = (C_o_L*L-HAB_o_L*L)^(-1)*(-DHABP_DC);    % auxiliary derivative \partial U_{t+1} / \partial C_t in steady state
@#endif
@#if UTIL == 1
DU_DC = (C_o_L*L-HAB_o_L*L)^(-ETAC);              % auxiliary derivative \partial U_t / \partial C_t in steady state
DUP_DC = (C_o_L*L-HAB_o_L*L)^(-ETAC)*(-DHABP_DC); % auxiliary derivative \partial U_{t+1} / \partial C_t in steady state
@#endif
LAM = XC^(-1)*(DU_DC + BETTA*DUP_DC); % foc C
%--------------------------------------------------------------------------
@#if MONPOL == 1
PIE = 1+PIEA/400;
pie = PIE;
R = PIE/BETTA;
r = R;
@#endif
%--------------------------------------------------------------------------
% Set remaining variables
y = Y_o_L*L;
c = C_o_L*L;
k = K_o_L*L;
iv = IV_o_L*L;
lam = LAM;
%--------------------------------------------------------------------------
end; % [steady_state_model] end
@#endif

% =========================================================================
% Initval
% =========================================================================
@#if LABOR == 1 && UTIL == 1
initval;
y = YBAR;
c = CBAR;
lam = LAMBAR;
iv = IVBAR;
k = KBAR;
rk = RKBAR;
q = QBAR;
a = ABAR;
@#if INVESTSHOCK == 1
upsilon = UPSILONBAR;
@#endif
@#if LABOR == 1
l = LBAR;
w = WBAR;
@#endif
@#if CAPUTIL == 1
uk = UKBAR;
@#endif
@#if MONPOL == 1
r = RBAR;
pie = PIEBAR;
@#endif
end; % [initval] end
@#endif

% =========================================================================
% Declare settings for shocks
% =========================================================================
shocks;
var epsa = 1;
@#if INVESTSHOCK == 1
var epsupsilon = 1;
@#endif
@#if MONPOL == 1
var epsr = 1;
@#endif
end; % [shocks] end


% =========================================================================
% Specify Priors
% =========================================================================
estimated_params;
% --------------------------------------------------------------------------------------------------------------------------------------------------
%PARAMETER_NAME, INITIAL_VALUE, LOWER_BOUND, UPPER_BOUND, PRIOR_SHAPE,   PRIOR_MEAN, PRIOR_STANDARD_ERROR, PRIOR_3RD_PARAMETER, PRIOR_4TH_PARAMETER;
% --------------------------------------------------------------------------------------------------------------------------------------------------
ALPHA,           0.3,           1e-8,        0.9999,      normal_pdf,    0.3,        0.05;
RA,              2,             1e-8,        10,          gamma_pdf,     2,          0.25;
DELTA,           0.025,         1e-8,        0.9999,      uniform_pdf,   ,           ,                     0,                   1;
RHOA,            0.5,           1e-8,        0.9999,      beta_pdf,      0.5,        0.1;
SIGA,            0.6,           1e-8,        10,          inv_gamma_pdf, 0.6,        2;
@#if MSADJ == 1
THETA,           1.5,           1e-8,        10,          gamma_pdf,     1.5,        0.75;
@#endif
@#if IAC > 0
KAPPA,           2,             1e-8,        10,          gamma_pdf,     2,          1.5;
@#endif
@#if CAPUTIL == 1
PSIK,            0.97,          1e-8,        0.9999,      uniform_pdf,   ,           ,                     0,                   1;
@#endif
@#if HABIT > 0
HBAR,            0.6,           1e-8,        1,           beta_pdf,      0.6,        0.1;
@#endif
@#if UTIL != 0
ETAC,            2,             1e-8,        10,          normal_pdf,    2,          0.3;
@#endif
@#if LABOR == 1
GAMMAL,          1.0052,        1e-8,        10,          uniform_pdf,   ,           ,                     0,                   10;
    @#if UTIL !=0
ETAL,            1.5,           1e-8,        10,          gamma_pdf,     1.5,        0.75;
    @#endif
@#endif
@#if INVESTSHOCK == 1
RHOUPSILON,      0.5,           1e-8,        0.9999,      beta_pdf,      0.5,         0.1;
SIGUPSILON,      0.6,           1e-8,        10,          inv_gamma_pdf, 0.6,         2;
@#endif
@#if MONPOL == 1
PIEA,            3.2,           0,           10,          gamma_pdf,     3.2,         2;
PSIPIE,          1.5,           1e-8,        10,          gamma_pdf,     1.5,         0.25;
PSIY,            0.125,         1e-8,        10,          gamma_pdf,     0.125,       0.1;
RHOR,            0.5,           1e-8,        0.9999,      beta_pdf,      0.5,         0.2;
SIGR,            0.2,           1e-8,        10,          inv_gamma_pdf, 0.2,         2;
@#endif
end; % [estimated_params] end


@#if COMPUTATIONS >= 0
% =========================================================================
% Steady-State, Checks and Diagnostics
% =========================================================================
steady(maxit=50000);                 % compute steady state given the starting values
resid;                  % check the residuals of model equations evaluated at steady state
check;                  % check Blanchard-Kahn conditions
model_diagnostics;      % check obvious model errors
disp(' ');
@#endif

@#if COMPUTATIONS == 1
% =========================================================================
% Lack of Identification Analysis
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
        identification(ar=10, no_identification_strength, parameter_set=calibration); % run identification analysis
        fprintf('**** Bruteforce VAROBS: %s ****' , strjoin(options_.varobs));
        identification(ar=10, no_identification_strength, parameter_set=calibration,checks_via_subsets=1, max_dim_subsets_groups=20); % run identification analysis
        fprintf('**** Monte Carlo Testing VAROBS: %s ****' , strjoin(options_.varobs));
        identification(ar=10, no_identification_strength, prior_mc=100, nograph, no_identification_minimal); % run identification analysis

        % save results into structure
        RESULTS{iset}{ii} = load([M_.fname, '/identification/', M_.fname, '_calibration_identif.mat'], 'ide_moments_point', 'ide_spectrum_point', 'ide_minimal_point');
        RESULTS{iset}{ii}.varobs = options_.varobs;
        
        % find out if KAPPA and/or THETA are among problematic parameter sets
        idx_KAPPA = find(estim_params_.param_vals(:,1)==find(contains(M_.param_names,'KAPPA')));
        idx_THETA = find(estim_params_.param_vals(:,1)==find(contains(M_.param_names,'THETA')));
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
                    if ~any(cellfun(@(x) nnz(ismember(x,[idx_KAPPA, idx_THETA])),ide.problpars(indx_sets), 'UniformOutput', 1))
                        TABLE{irun,jide} = '$\checkmark$'; % both theta and kappa are identified, but other parameters are not
                    else
                        TABLE{irun,jide} = '$[';
                        if any(cellfun(@(x) nnz(ismember(x,idx_KAPPA)),ide.problpars(indx_sets), 'UniformOutput', 1))
                            TABLE{irun,jide} = [TABLE{irun,jide} '\kappa ']; % kappa is within unidentified sets
                        end
                        if any(cellfun(@(x) nnz(ismember(x,idx_THETA)),ide.problpars(indx_sets), 'UniformOutput', 1))
                            TABLE{irun,jide} = [TABLE{irun,jide} '\theta ']; % theta is within unidentified sets
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
    
% save LaTex code as file and compile it
filenam = sprintf('table.tex');
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
save('../simdat.mat','y','c','iv','rk','k','lam','q','a');
@#if CAPUTIL == 1
save('../simdat.mat','uk','-append');
@#endif
@#if LABOR == 1
save('../simdat.mat','l','w','-append');
@#endif
@#if MONPOL == 1
save('../simdat.mat','r','pie','-append');
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
    elseif jj==length(optimtypevals)
        eval(['options_.mode_compute=' num2str(numopt) ';']) ;        
        options_.gmhmaxlik.target=0.2;
        options_.gmhmaxlik.iterations = 1;  % [default=3] to call repeatedly the new optimization routine. Improves the estimates of the posterior covariance matrix and of the posterior mode.
        estimation(optim=('AcceptanceRateTarget',0.2),mode_check,mode_file=newmode,cova_compute=1,mh_replic=0,datafile=simdat,first_obs=@{FIRSTOBS},nobs=@{SAMPLESIZE},tex);
    end
end
movefile([M_.fname '_mode.mat'], 'finalmode.mat'); % rename mode file
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
