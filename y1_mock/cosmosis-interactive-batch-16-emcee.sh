#!/bin/bash

# N.B.: Our code expects that each MPI rank can see all the GPUs on each node.
# The code handles making sure that each MPI rank talks only to a single GPU.
# We must use the --gpus-per-task=1 and --gpu-bind=none flags for SLURM because
# of this.
#
# Make sure to keep 'ntasks' = 4 * 'nodes', because always want 4 ranks per node.
# Use -q debug only for 4 or fewer nodes
# Use -q regular for more than 4 nodes, up to 128 nodes

#SBATCH -A des_g
#SBATCH -C gpu
#SBATCH -G 16
#SBATCH -q debug
#SBATCH -t 00:30:00
#SBATCH --nodes=4
#SBATCH --ntasks=16
#SBATCH --ntasks-per-node=4
#SBATCH -c 32
#SBATCH --gpus-per-task=1
#SBATCH --gpu-bind=none
#SBATCH --mail-type=begin,end,fail
#SBATCH --mail-user=jesteves@umich.edu

# working directory, use high performance shared area on perlmutter

export TOP_DIR=/global/common/software/des/jesteves
export COSMOSIS_REPO_DIR=${TOP_DIR}/cosmosis
export CSL_DIR=${TOP_DIR}/cosmosis-standard-library

export INTEGRATION_TOOLS_DIR=${TOP_DIR}/y3_pipe_under
export Y3PIPE_DIR=${TOP_DIR}/y3_cluster_cpp
export Y3_CLUSTER_WORK_DIR=${Y3PIPE_DIR}/release-build
export Y3_CLUSTER_CPP_DIR=${Y3PIPE_DIR}
export COSMOSIS_STANDARD_LIBRARY=${CSL_DIR}

export CUBA_DIR=${INTEGRATION_TOOLS_DIR}/cuba
export CUBA_CPP_DIR=$INTEGRATION_TOOLS_DIR/cubacpp
export GPU_INT_DIR=${INTEGRATION_TOOLS_DIR}/gpuintegration 

# We set OMP_NUM_THREADS to avoid oversubscribing the CPU cores.
# This value has not been carefully tuned.
export OMP_NUM_THREADS=4

# Now set up CosmoSIS
source ${COSMOSIS_REPO_DIR}/setup-cosmosis-nersc /global/common/software/des/common/Conda_Envs/y3cl_je

# And let's work in our own individual y3_cluster_cpp directory
export MY_TOP_DIR=/global/common/software/des/$(id -un)
export Y3PIPE_DIR=${MY_TOP_DIR}/y3_cluster_cpp
export Y3_CLUSTER_CPP_DIR=${Y3PIPE_DIR}
export Y3_CLUSTER_WORK_DIR=${Y3PIPE_DIR}/release-build
cd ${Y3PIPE_DIR}/y1_mock/
export SLURM_CPU_BIND="cores"

srun -n ${SLURM_NTASKS} cosmosis --mpi mock.ini
