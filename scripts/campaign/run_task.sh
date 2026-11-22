#!/bin/bash
# Worker wrapper for the gate-0 array jobs.
#   run_task.sh <python script> [args...]
#
# MPI4PY_RC_INITIALIZE=0 is the important line.  DRAC's netCDF4 wheel is linked
# against MPI and pulls in mpi4py, which calls MPI_Init on import.  Inside a
# single-task Slurm step that aborts with "An error occurred in MPI_Init_thread
# on a NULL communicator".  These are serial jobs that only need netCDF for
# serial I/O, so we tell mpi4py not to initialise MPI at all.
set -euo pipefail
source "$HOME/envs/meromixis-env/activate-meromixis"
export MPI4PY_RC_INITIALIZE=0
export MPI4PY_RC_FINALIZE=0
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
exec python "$@"
