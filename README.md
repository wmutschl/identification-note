[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3192196.svg)](https://doi.org/10.5281/zenodo.3192196)

# Replication files for "The effect of observables, functional specifications, model features and shocks on identification in linearized DSGE models"
- Authors: Sergey Ivashchenko and Willi Mutschler (corresponding author)
- Published in [Economic Modelling](https://doi.org/10.1016/j.econmod.2019.09.039)

### General notes
- <u>**For all cases you can simply check the Dynare log files and generated pdfs instead of (time-costly) rerunning all models.**</u>
- Note that you need at least Dynare 4.6-unstable. All results were generated with the Version from May 6, 2019 (commit 0375dbe29b40ac1f70963ccfd8e429ef99fea774 on [Dynare Master Branch](https://git.dynare.org/Dynare/dynare). You can find this specific version in the folder `utils`.
- If you spot mistakes please open an issue or write an email to willi@mutschler.eu

### Updates
- Small bugfixes
- Rerun all identification rank checks for prior mean and 100 draws from prior domain
- Added measurement errors to AnScho model
- Rerun all identification strength scenarios and updated tables

### Folder `utils`
- `dynare-4.6-unstable-0375dbe.tar.xz` and `dynare-4.6-unstable-0375dbe-win.exe` source files and Windows Setup file of the Dynare version used for replication. Any newer Dynare 4.6-unstable version should work as well. There might be slight differences due to random numbers and bugfixes. We would always advise to use the latest dynare 4.6 version from [https://dynare.org](https://dynare.org)
- `AnSchoModTheBuilder.mod`: Dynare mod file for all variants of the monetary model. This mod file is copied into all folders.
- `KimModTheBuilder.mod`: Dynare mod file for all variants of the investment adjustment costs model. This mod file is copied into all folders.
- `LatexTable.m` function that converts Matlab tables to Latex tables. Original authors are Eli Duenisch and Pascal E. Fortin.
- `SetupForParallel.ini` configuration file for running Dynare in parallel on a local computer with eight cores.
- `SetupForParallelCluster.ini` configuration file for running Dynare in parallel on the PALMA cluster of the University of Münster with 72 cores.

### Folder `Kim/Lack_of_Identification`
- Contains replication files and results for the theoretical identification analysis by analyzing the rank criteria of *Iskrev (2010)*, *Komunjer and Ng (2011)* and *Qu and Tkachenko (2012)*
- We consider the following model scenarios
    - BASELINE: The original model of *Kim (2003, JEDC)*
    - CAPUTILIZATION: Baseline with capital utilization
    - CRRA_EXTERNALHABIT: Baseline with CRRA utility function and external consumption habit formation
    - CRRA_INTERNALHABIT: Baseline with CRRA utility function and internal consumption habit formation
    - CRRA_NOHABIT: Baseline with CRRA utility function
    - EXTERNALHABIT: Baseline with external consumption habit formation
    - INTERNALHABIT: Baseline with internal consumption habit formation
    - INVESTSPECSHOCK: Baseline with investment-specific technological shock
    - LABOR: Baseline with leisure/labor choice
    - MONPOL: Baseline with monetary policy
- Model folders
    - The first part of the name of the model folder indicates the scenario, whereas the last part indicates the specification of intertemporal investment adjustment costs (*level* or *growth*)
    - `RUN_DYNARE.m` function that sets the Dynare macros and needs to be run in each folder
    - `KimModTheBuilder.log` contains the log file of the identification analysis (here you can simply check the run instead of time consuming rerunning everything)
    - `table.pdf` contains a table with a summary of the results
- `Replicate_Kim_Lack_of_Identification.m`
    - Recreates/changes a specific or all model variants. You then have to go into each model folder and run `RUN_DYNARE.m`. This will recreate all log files as well as pdf files. Be careful as this function deletes all previous results and makes the folders empty.

### Folder `Kim/Strength_of_Identification`
- Contains replication files and results for the Bayesian learning rate indicator of Koop, Pesaran, and Smith (2013). To this end, we need to simulate data and then estimate the posterior precision for different (growing) sample sizes.
- We consider the following model scenarios with either the *growth* or *level* specification of intertemporal investment adjustment costs
    - Baseline: The original model of *Kim (2003, JEDC)*
    - Investshock: Baseline with investment-specific technological shock
- `RunKimStrength.m` creates all model folders and contains all details for the Dynare macros needed for all model variants, the other matlab files starting with `Replicate...` simply set different sample sizes.
- If you want to rerun an estimation for a certain model variant with a certain sample size, e.g. for the Baseline growth scenario and 100 observations simply run `ReplicateKimBaselineGrowth_100.m`. This will create a model folder (`BaselineGrowth/100`) with the following files:
    - `currentmodelcall.txt` the exact invocation of Dynare used in this folder
    - `KimModTheBuilder.log` is the log file of the estimation (simply check this instead of rerunning everything)
    - `KimModTheBuilder_TeX_binder.pdf` is a summary of all graphs (in particular convergence diagnostics) produced by Dynare during the estimation
    - The average posterior precisions are saved in MAT files in the parent folder (e.g. `weakresults_mcmc_100.mat` in `BaselineGrowth` folder)
- After all estimations are done, `Make_Kim_Latex_Tables_Strength.m` creates Latex tables in each Scenario folder into a folder labled `tables` from the `weakresults_mcmc_...` mat files.

### Folder `AnScho/Lack_of_Identification`
- Contains replication files and results for the theoretical identification analysis by analyzing the rank criteria of *Iskrev (2010)*, *Komunjer and Ng (2011)* and *Qu and Tkachenko (2012)*
- We consider the following model scenarios
    - BASELINE: The original model of *An and Schorfheide (2007, Econometric Reviews)*
    - INDEXATION: Baseline with partial inflation indexation
    - PREFSHOCK: Baseline with a preference shock on the discount factor shifter
    - INDEXATION and PREFSHOCK: Baseline with both partial inflation indexation and preference shock
- Model folders
    - The first part of the name of the model folder indicates the scenario, whereas the last part indicates the specification of monetary policy rule (*FLEX*, *GROWTH*, *STEADYSTATE* or *SW*)
    - We also added a version with measurement errors on YGR, INFL, and INT, these folders end with *MEASERR*
    - `RUN_DYNARE.m` function that sets the Dynare macros and needs to be run in each folder
    - `AnSchoTheBuilder.log` contains the log file of the identification analysis (simply check this instead of rerunning the estimation)
    - `table.pdf` contains a table with a summary of the results
- `Replicate_AnScho_Lack_of_Identification.m`
    - Recreates/changes a specific or all model variants. You then have to go into each model folder and run `RUN_DYNARE.m`. This will recreate all log files as well as pdf files. Be careful as this function deletes all previous results and makes the folders empty.

### Folder `AnScho/Strength_of_Identification`
- Contains replication files and results for the Bayesian learning rate indicator of Koop, Pesaran, and Smith (2013). We need to simulate data and then estimate the posterior precision for different (growing) sample sizes
- We consider the following model scenarios:
    - `Baseline`: Original model with flex-price monetary policy rule
    - `Indexation`: Original model with partial inflation indexation and flex-price monetary policy rule
    - `MonPolSteadyStateGap`: Original model with steady-state monetary policy rule
    - `Prefshock`: Original model with preference shock and flex-price monetary policy rule
- `RunAnSchoStrength.m` creates all model folders and contains all details for the Dynare macros needed for all model variants, the other matlab files starting with `Replicate...` simply set different sample sizes.
- If you want to rerun an estimation for a certain model variant with a certain sample size, e.g. for the Baseline scenario and 100 observations simply run `ReplicateASBaseline_100.m`. This will create a model folder (`Baseline/100`) with the following files:
    - `currentmodelcall.txt` the exact invocation of Dynare used in this folder
    - `AnSchoModTheBuilder.log` is the log file of the estimation (simply check this instead of rerunning the estimation)
    - `AnSchoModTheBuilder_TeX_binder.pdf` is a summary of all graphs (in particular convergence diagnostics) produced by Dynare during the estimation
    -  The average posterior precisions are saved in MAT files in the parent folder (e.g. `weakresults_mcmc_100.mat` in `Baseline` folder)
- After all estimations are done, `Make_AnScho_Latex_Tables_Strength.m` creates Latex tables in each Scenario folder into a folder labled `tables` from the `weakresults_mcmc_...` mat files.
