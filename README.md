# Numerical method for approximately optimal solutions of two-stage distributionally robust optimization with marginal constraints

+ By Ariel Neufeld and Qikun Xiang
+ Article link (arXiv): https://arxiv.org/abs/2205.05315

# Description of files

+ func/mixtruncnorm/      contains functions related to mixture of truncated normal distributions
    - norm\_partialexp.m: function used for computing a specific form of expectation related to a normal distribution
    - mixtruncnorm\_partialexp.m: function used for computing a specific form of expectation related to a mixture of truncated normal distributions
    - mixtruncnorm\_invcdf.m: function used for computing the inverse cumulative distribution function of a mixture of truncated normal distributions via the bisection method
    - mixtruncnorm\_momentset.m: function used for computing the integrals of the functions in an interpolation function basis (which characterizes a moment set) with respect to a mixture of truncated normal distributions
    - mixtruncnorm\_momentset\_construct.m: function used for iteratively constructing a moment set surrounding a mixture of truncated normal distributions with a given number of knots

+ func/momentset1d/      contains functions related to moment set surrounding one-dimensional distribution with bounded support
    - momentset1d\_basisfunc\_bounded.m: function used for computing the values of the functions in an interpolation function basis for given inputs
    - momentset1d\_measinit\_bounded.m: function that returns a discrete measure from a moment set characterized by a given interpolation function basis

+ func/CPWA/       contains functions related to univariate continuous piece-wise affine (CPWA) functions
    - CPWA\_1d\_conjugate.m: function used for computing the convex conjugate of a one-dimensional convex piecewise affine functions with compact domain
    - CPWA\_1d\_convexenv.m: function used for computing the convex envelope function of a one-dimensional piece-wise affine function with compact domain
    - CPWA\_sep\_concavemin\_MILP\_gurobi.m: function that returns a Gurobi model which formulates the problem of minimizating a separable concave function over a polytope into a mixed-integer linear programming (MILP) problem

+ func/reassembly/       contains functions related to coupling and reassembly
    - comonotone\_coupling.m: function used for computing the comonotone coupling of discrete measures (i.e., a joint distribution formed with the given set of discrete marginals and the comonotone copula)
    - discretecopula\_samp.m: function used for generating independent samples from the copula of a given discrete measure
    - reassembly\_increment.m: function used for computing the reassembly of a given discrete measure with a given set of discrete marginals

+ func/DRO/       contains functions used in the algorithm for approximately solving two-stage distributionally robust optimization (DRO) problems with marginal constraints
    - DRO\_cutplane.m: implementation of the cutting plane algorithm for solving linear semi-infinite programming formulations of relaxed two-stage DRO problems
    - DRO\_feascons.m: function that returns constraints (also known as feasibility cuts) in the cutting plane algorithm
    - DRO\_primal\_approx.m: function used for approximately computing a lower bound for the optimal value of the two-stage DRO problem with marginal constraints via Monte Carlo integration
    - DRO\_primal\_tighten\_gurobi.m: function that returns a Gurobi model which corresponds to a linear programming problem used for generating/improving a lower bound for the optimal value of the relaxed two-stage DRO problem (see Supplementary.pdf)
    - DRO\_stage2\_gurobi.m: function that returns a Gurobi model which formulates the second-stage problem in DRO for a given first-stage decision and a given vector of uncertain quantities
    - SP\_approx.m: function that formulates and solves the sample average approximation (SAA) of a stochastic programming problem given a collection of samples (only used for checking the feasibility of the DRO problem)

+ func/DRO/problems			contains functions used for formulating concrete two-stage DRO problems
    - DRO\_task\_scheduling.m: function that returns the DRO problem formulation of a task scheduling problem
    - DRO\_multi\_product\_assembly.m: function that returns the DRO problem formulation of a multi-product assembly problem
    - DRO\_supply\_chain\_network.m: function that returns the DRO problem formulation of a supply chain network design problem with edge failure
    - DRO\_supply\_chain\_network\_marginv.m: utility function for computing the inverse cdf transform in the supply chain network problem with edge failure
    - SP\_supply\_chain\_network.m: function that returns the stochastic programming problem formulation of a supply chain network design problem with edge failure (only used for checking the feasibilty of the corresponding DRO problem)

+ exp/            contains the scripts to run the numerical experiments (see below for detailed instructions)

+ utils/          contains external libraries
    - utils/tight\_subplot/:             used for creating figures with narrow margins
    - utils/parfor\_progress/:           used for printing a progress bar in the console

# Instructions to run the numerical experiments

## Configurations

+ All folders and subfolders must be added to the MATLAB search path. 
+ Gurobi optimization (version 9.5.0 or above) must be installed on the machine and relevant files must be added to the search path. 

## The task scheduling problem

### Step 1: generate the input files
+ Run exp/DRO\_scheduling/DRO\_scheduling\_exp\_prepare1.m to generate a file exp/DRO\_scheduling/exp\_inputs.mat containing the two-stage DRO specifications and the approximation scheme.
+ Run exp/DRO\_scheduling/DRO\_scheduling\_exp\_prepare2.m to generate a file exp/DRO\_scheduling/exp\_invcdf.mat containing the necessary inputs for approximating the inverse cumulative distribution functions of the marginals.

### Step 2: compute the upper and lower bounds for the optimal value of the two-stage DRO problem with marginal constraints
+ Run exp/DRO\_scheduling/DRO\_scheduling\_exp\_run1\_UB.m to compute the upper bounds by solving linear semi-infinite programming problems via the cutting plane algorithm. This will create an output file exp/DRO\_scheduling/exp\_rst\_UB.mat.
+ Run exp/DRO\_scheduling/DRO\_scheduling\_exp\_run2\_LB.m to approximately compute the lower bounds via Monte Carlo integration with respect to a constructed reassembly. This will create an output file exp/DRO\_scheduling/exp\_rst\_LB.mat.

### Step 3: plot the results
+ Run exp/DRO\_scheduling/DRO\_scheduling\_exp\_plot\_results.m to plot the upper and lower bounds, the comparison of the differences between the bounds and their theoretical upper bounds, and the approximately optimal scheduled duration of the tasks.


## The multi-product assembly problem
### Step 1: generate the input files
+ Run exp/DRO\_assembly/DRO\_assembly\_exp\_prepare1.m to generate a file exp/DRO\_assembly/exp\_inputs.mat containing the two-stage DRO specifications and the approximation scheme.
+ Run exp/DRO\_assembly/DRO\_assembly\_exp\_prepare2.m to generate a file exp/DRO\_assembly/exp\_invcdf.mat containing the necessary inputs for approximating the inverse cumulative distribution functions of the marginals.

### Step 2: compute the upper and lower bounds for the optimal value of the two-stage DRO problem with marginal constraints
+ Run exp/DRO\_assembly/DRO\_assembly\_exp\_run1\_UB.m to compute the upper bounds by solving linear semi-infinite programming problems via the cutting plane algorithm. This will create an output file exp/DRO\_assembly/exp\_rst\_UB.mat.
+ Run exp/DRO\_assembly/DRO\_assembly\_exp\_run2\_LB.m to approximately compute the lower bounds via Monte Carlo integration with respect to a constructed reassembly. This will create an output file exp/DRO\_assembly/exp\_rst\_LB.mat.

### Step 3: plot the results
+ Run exp/DRO\_assembly/DRO\_assembly\_exp\_plot\_results.m to plot the upper and lower bounds.

## The supply chain network design problem with edge failure
### Step 1: generate the input files
+ Run exp/DRO\_network/DRO\_network\_exp\_prepare1.m to generate a file exp/DRO\_network/exp\_inputs.mat containing the two-stage DRO specifications and the approximation scheme.
+ Run exp/DRO\_network/DRO\_network\_exp\_prepare2.m to generate a file exp/DRO\_network/exp\_invcdf.mat containing the necessary inputs for approximating the inverse cumulative distribution functions of the demand distributions.

### Step 2: compute the upper and lower bounds for the optimal value of the two-stage DRO problem with marginal constraints
+ Run exp/DRO\_network/DRO\_network\_exp\_run1\_UB.m to compute the upper bounds by solving linear semi-infinite programming problems via the cutting plane algorithm. This will create an output file exp/DRO\_network/exp\_rst\_UB.mat.
+ Run exp/DRO\_network/DRO\_network\_exp\_run2\_LB.m to approximately compute the lower bounds via Monte Carlo integration with respect to a constructed reassembly. This will create an output file exp/DRO\_network/exp\_rst\_LB.mat.

### Step 3: plot the results
+ Run exp/DRO\_network/DRO\_network\_exp\_plot\_network.m to plot network structure along with the approximately optimal investments for the processing facilities.
+ Run exp/DRO\_network/DRO\_network\_exp\_plot\_results.m to plot the upper and lower bounds.


# About the supplementary material
The document Supplementary.pdf contains explanations about how the initial set of constraints is generated for the cutting-plane algorithm when solving the linear semi-infinite programming problem.
