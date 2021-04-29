# TAMMBER
### Temperature Accelerated Markov Model construction with Bayesian Estimation of Rates
TAMMBER is a heavily modified variant of the [ParSplice](https://gitlab.com/exaalt/parsplice.git) code (see below).

Designed for massively parallel deployment, TAMMBER optimally manages hundreds to thousands of workerso 
performing accelerated molecular dynamics (TAD) and minimum energy path routines (NEB). 
This sampling data is processed using Bayesian techniques to yeild uncertainty-controlled kMC/Markov models
of complex atomistic target systems, with minimal end-user involvement. [Recent](https://www.nature.com/articles/s41524-020-00463-8) 
developements allow for the autonomous construction and convergence of arbitrarily complex diffusion tensors.

To be efficient at the peta- or exa-scale requires optimal worker management. TAMMBER achieves this 
by treating as-yet-unseen configuration space as an absorbing sink, then calculating (given the current information)
where additional sampling effort will maximally increase the time-to-absorbtion through discovery or reduction of uncertainty.

TAMMBER can efficiently scale to 100,000+ cores, though can be profitibly employed on 200+. See the publications below for more detail.

Questions / bugs ? Raise a github issue or email swinburne "at" cinam "dot" univ-mrs "dot" fr

##  User Guide
### - [Installation Instructions](INSTALL.md)
### - [Getting Started](EXAMPLE.md)
   - *If you are using TAMMBER for diffusion studies please read the "Cluster Definitions For Diffusion" section!*
   - *It is highly recommended to run the `tammber-md` test routine first!*
### - [Output Analysis](process/Diffusion_Model_Example.ipynb) :  Run online with Binder: [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/tomswinburne/tammber/HEAD?filepath=process%2FDiffusion_Model_Example.ipynb)
   - *If the Binder server is slow the notebook can be viewed [here](process/Diffusion_Model_Example.ipynb), or run the notebook locally after downloading this repository.*

## Coming Soon
  - Updated analysis scripts with custom point groups and disconnectivity graphs
  - More detailed analysis of NEB simulations "on-the-fly" (specifically for multiple minima)

## Publications:
[1] Thomas D Swinburne and Danny Perez, *Self-optimized construction of transition rate matrices from accelerated atomistic simulations with Bayesian uncertainty quantification*, Physical Review Materials 2018, [preprint](https://arxiv.org/abs/1803.05273)
   *Original publication detailing the optimal worker management and Bayesian uncertainty quantification*
[2] Thomas D Swinburne and Danny Perez, *Automated Calculation of Defect Transport Tensors*, NPJ Computational Materials, 2020. [article](https://www.nature.com/articles/s41524-020-00463-8)
   *Incorporation of crystal symmetries to efficiently and autonomously calculate diffusion tensors. See the [Getting Started](EXAMPLE.md) section*




--------------------------------------------------------------------------------

# ParSplice code
The [ParSplice](https://gitlab.com/exaalt/parsplice.git) code implements the Parallel Trajectory Splicing algorithm described in [3]. This method is part of the Accelerated Molecular Dynamics family of techniques developed in Los Alamos National Laboratory over the last 16 years. These methods aim at generating high-quality trajectories of ensembles of atoms in materials. ParSplice uses multiple independent replicas of the system in order to parallelize the generation of such trajectories in the time domain, enabling simulations of systems of modest size over very long timescales. ParSplice includes capabilities to store configurations of the system, to generate and distribute tasks across a large number of processors, and to harvest the results of these tasks to generate long trajectories. ParSplice is a management layer that orchestrate large number of calculations, but it does not perform the actual molecular dynamics itself; this is done by external molecular dynamics engines.

[3] Danny Perez, Ekin D Cubuk, Amos Waterland, Efthimios Kaxiras, Arthur F Voter, *Long-time dynamics through parallel trajectory splicing*, Journal of chemical theory and computation 12, 18 (2015)
