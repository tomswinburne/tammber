#!/bin/bash
#SBATCH -J job_name
#SBATCH --nodes=10
#SBATCH --ntasks-per-node=28
#SBATCH --time=0:15:00
#SBATCH --exclusive
#SBATCH --mem=10GB # 
#SBATCH --output output.slurm 
set -e

# Running TAMMBER on 28*10=280 nodes for 15 minutes as a test run

# Load required modules (boost required if dynamic linking on some architectures)
module load intel/19.4 intelmpi/2019.4.243  boost/1_69_0
ulimit -s unlimited

TAMMBER_WORK_MANAGERS=1 # need one manager for every ~500 workers
TAMMBER_CORES_PER_WORKER=1 # one core per ~20k atoms with typical LAMMPS scaling
TAMMBER_PAD_NODES=0 # only required for very large jobs
TAMMBER_SEED=2345 # random number seed
srun -l /path/to/tammber/repo/build/tammber ${SLURM_JOB_NUM_NODES} ${SLURM_CPUS_ON_NODE} ${TAMMBER_CORES_PER_WORKER} ${TAMMBER_WORK_MANAGERS} ${TAMMBER_PAD_NODES} ${TAMMBER_SEED} 


