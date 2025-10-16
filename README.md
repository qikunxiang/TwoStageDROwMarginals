# Numerical method for approximately optimal solutions of two-stage distributionally robust optimization with marginal constraints

+ By Ariel Neufeld and Qikun Xiang
+ Article link (arXiv): https://arxiv.org/abs/2205.05315

## Description of files

+ **instance/** contains classes for describing instances of the two-stage distributionally robust optimization (DRO) problems with marginal constraints
	- **TSDROMInstance.m**: abstract class for two-stage DRO problems with marginal constraints
	- **TSDROMUncertain.m**: auxiliary class for handling the marginal constraints in a two-stage DRO problems and related computations
	- **TSDROMScheduling.m**: class for instances of the distributionally robust task scheduling problem
	- **TSDROMMultiProdAssem.m**: class for instances of the distributionally robust multi-product assembly problem
	- **TSDROMSupplyNet.m**: class for instances of the distributionally robust supply chain network design problem with edge failure
	- **TSDROMSupplyNetUncertain.m**: auxiliary class for handling the marginal constraints in the distributionally robust supply chain network design problem with edge failure; it is used to facilitate the initialization of the discrete and continuous marginals

+ **probmeas/** contains classes related to probability measures
	- **TSDROMProb.m**: abstract class for probability measures supported in a one-dimensional compact interval
	- **TSDROMProbBernoulli.m**: class for probability laws of Bernoulli distributions, i.e., probability measures consisting of two Dirac measures at 0 and 1
	- **TSDROMProbPWAwAtom.m**: class for probability measures with (possibly discontinuous) piece-wise affine (PWA) density functions mixed with discrete measures
	- **TSDROMProbTruncMixNorm.m**: class for probability measures with truncated mixture of normal density functions

+ **solver/** contains classes for the components of the solver
	- **LSIPMinCuttingPlaneAlgo.m**: abstract class for cutting-plane algorithms used for solving linear semi-infinite programming (LSIP) minimization problems
    - **TSDROMDualLSIPSolver.m**: class for the concrete implementation of the cutting-plane algorithm for solving the dual LSIP problem
	- **TSDROMNewVariableSolver.m**: abstract class for modeling the global maximization problem in order to identify new points in the support of the dual probability measure
	- **TSDROMNewVariableSolverExtrapolation.m**: abstract class for the global maximization problem where the boundaries of extrapolation for each Kantorovich potential function need to be handled iteratively
	- **TSDROMNewVariableINCSolver.m**: class for the concrete implementation of the global maximization problem via formulating the problem into a mixed-integer programming (MIP) problem with the incremental model
	- **TSDROMNewVariableLOGSolver.m**: class for the concrete implementation of the global maximization problem via formulating the problem into an MIP problem with the logarithmic model
	- **TSDROMNewVariablePWLSolver.m**: class for the concrete implementation of the global maximization problem using Gurobi's piece-wise linear objective function utilities
	- **TSDROMIterativeSolver.m**: class for the concrete implementation of the iterative numerical algorithm for computing approximately optimal decisions as well as approximate worst-case probability measures
	- **TSDROMWorstCaseSampler.m**: auxiliary class for efficiently generating independent random samples from a computed approximate worst-case probability measure


+ **experiments/** contains the scripts to run the numerical experiments (see below for detailed instructions)

+ **utils/** contains external libraries
   - **tight\_subplot/**: used for creating figures with narrow margins

## Instructions to run the numerical experiments

### Configurations

+ All folders and subfolders must be added to the MATLAB search path. 
+ Gurobi optimization (version 12.0.3 or above) must be installed on the machine and the relevant files must be added to the search path. 
+ In the file **experiments/global\_config.m**, the values of *CONFIG.SAVEPATH_ROOT* and *CONFIG.LOGPATH_ROOT* can be modified to change the directories where save files and log files are saved. 

### Experiment: task scheduling

#### Step 1: generate the input file
+ Run **experiments/Scheduling\_Exp/Scheduling\_exp\_prepare.m** to generate a file containing the settings of the distributionally robust task scheduling problem.

#### Step 2: compute an approximately optimal decision and an approximate worst-case probability measure
+ Run **experiments/Scheduling\_Exp/Scheduling\_exp\_run.m** to execute the numerical algorithm.
+ Run **experiments/Scheduling\_Exp/Scheduling\_exp\_print\_results.m** to print the outputs of the numerical algorithm and some statistics related to the execution of the algorithm. 
+ Run **experiments/Scheduling\_Exp/Scheduling\_exp\_examine\_measure.m** to examine the theorized property of the computed approximate worst-case probability measure. 

#### Step 3: visualize the results
+ Run **experiments/Scheduling\_Exp/Scheduling\_exp\_plot\_sparsity.m** to plot the scatter plots of the computed upper and lower bounds and the computed sub-optimality estimates against the cardinalities of the supports of the dual measures.
+ Run **experiments/Scheduling\_Exp/Scheduling\_exp\_plot\_decision.m** to plot the computed approximately optimal scheduled durations of the tasks. 
+ Run **experiments/Scheduling\_Exp/Scheduling\_exp\_plot\_potential.m** to plot four computed subdifferential mappings.
+ Run **experiments/Scheduling\_Exp/Scheduling\_exp\_plot\_dependence.m** to plot the pair-wise scatter plots of four marginals of the computed approximate worst-case probability measure.


### Experiment: multi-product assembly

#### Step 1: generate the input file
+ Run **experiments/MultiProdAssem\_Exp/MultiProdAssem\_exp\_prepare.m** to generate a file containing the settings of the distributionally robust multi-product assembly problem.

#### Step 2: compute an approximately optimal decision and an approximate worst-case probability measure
+ Run **experiments/MultiProdAssem\_Exp/MultiProdAssem\_exp\_run.m** to execute the numerical algorithm.
+ Run **experiments/MultiProdAssem\_Exp/MultiProdAssem\_exp\_print\_results.m** to print the outputs of the numerical algorithm and some statistics related to the execution of the algorithm. 

#### Step 3: visualize the results
+ Run **experiments/MultiProdAssem\_Exp/MultiProdAssem\_exp\_plot\_sparsity.m** to plot the scatter plots of the computed upper and lower bounds and the computed sub-optimality estimates against the cardinalities of the supports of the dual measures.
+ Run **experiments/MultiProdAssem\_Exp/MultiProdAssem\_exp\_plot\_potential.m** to plot four computed subdifferential mappings.
+ Run **experiments/MultiProdAssem\_Exp/MultiProdAssem\_exp\_plot\_dependence.m** to plot the pair-wise scatter plots of four marginals of the computed approximate worst-case probability measure.


### Experiment: supply chain network design with edge failure

#### Step 1: generate the input file
+ Run **experiments/SupplyNet\_Exp/SupplyNet\_exp\_prepare.m** to generate a file containing the settings of the distributionally robust supply chain network design problem with edge failure.

#### Step 2: compute an approximately optimal decision and an approximate worst-case probability measure
+ Run **experiments/SupplyNet\_Exp/SupplyNet\_exp\_run.m** to execute the numerical algorithm.
+ Run **experiments/SupplyNet\_Exp/SupplyNet\_exp\_print\_results.m** to print the outputs of the numerical algorithm and some statistics related to the execution of the algorithm. 

#### Step 3: visualize the results
+ Run **experiments/SupplyNet\_Exp/SupplyNet\_exp\_plot\_sparsity.m** to plot the scatter plots of the computed upper and lower bounds and the computed sub-optimality estimates against the cardinalities of the supports of the dual measures.
+ Run **experiments/SupplyNet\_Exp/SupplyNet\_exp\_plot\_network.m** to plot the supply chain network configuration along with the computed approximately optimal processing capabilities.
+ Run **experiments/SupplyNet\_Exp/SupplyNet\_exp\_plot\_potential.m** to plot four computed subdifferential mappings.
+ Run **experiments/SupplyNet\_Exp/SupplyNet\_exp\_plot\_dependence.m** to plot the pair-wise scatter plots of four marginals of the computed approximate worst-case probability measure.