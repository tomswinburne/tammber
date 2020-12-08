# TAMMBER 
### Temperature Accelerated Markov Model construction with Bayesian Estimation of Rates

##  [Installation Instructions](INSTALL.md)

##  [Getting Started](EXAMPLE.md)

## [Output Analysis](process/Diffusion_Model_Example.ipynb) :  Run online with Binder: [![Launch](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gl/tomswinburne%2Fparsplice/master?filepath=process%2FDiffusion_Model_Example.ipynb)
*If the Binder server is slow the notebook can be viewed [here](process/Diffusion_Model_Example.ipynb), or run the notebook locally after downloading this repository.*

Questions? swinburne "at" cinam "dot" univ-mrs "dot" fr
--------------------------------------------------------------------------------

TAMMBER is a heavily modified variant of the [ParSplice](https://gitlab.com/exaalt/parsplice.git) code (see below). Whilst ParSplice constructs a single long trajectory by splicing many short MD segments, TAMMBER uses TAD-MD and barrier calculations to build rate matrices with a novel measure of the "missing" rate from seen to unseen regions of configuration space.

The work management is controlled through an absorbing Markov chain; the goal is to maximize the time to absorption from a given starting distribution, with statewise absorption rates calculated through Bayesian analysis of the trajectory segments

[1] Thomas D Swinburne and Danny Perez, *Self-optimized construction of transition rate matrices from accelerated atomistic simulations with Bayesian uncertainty quantification*, Physical Review Materials 2018, [preprint](https://arxiv.org/abs/1803.05273)


New code incorporates crystal symmetries using graph isomorphisms to efficiently and autonomously calculate diffusion matricies 

[2] Thomas D Swinburne and Danny Perez, *Automated Calculation of Defect Transport Tensors*, Accepted in NPJ Computational Materials, 2020. [preprint](https://arxiv.org/abs/2003.07752)

--------------------------------------------------------------------------------

# ParSplice code
The [ParSplice](https://gitlab.com/exaalt/parsplice.git) code implements the Parallel Trajectory Splicing algorithm described in [3]. This method is part of the Accelerated Molecular Dynamics family of techniques developed in Los Alamos National Laboratory over the last 16 years. These methods aim at generating high-quality trajectories of ensembles of atoms in materials. ParSplice uses multiple independent replicas of the system in order to parallelize the generation of such trajectories in the time domain, enabling simulations of systems of modest size over very long timescales. ParSplice includes capabilities to store configurations of the system, to generate and distribute tasks across a large number of processors, and to harvest the results of these tasks to generate long trajectories. ParSplice is a management layer that orchestrate large number of calculations, but it does not perform the actual molecular dynamics itself; this is done by external molecular dynamics engines.

[3] Danny Perez, Ekin D Cubuk, Amos Waterland, Efthimios Kaxiras, Arthur F Voter, *Long-time dynamics through parallel trajectory splicing*, Journal of chemical theory and computation 12, 18 (2015)
