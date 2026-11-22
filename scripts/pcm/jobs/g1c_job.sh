#!/bin/bash
#SBATCH --account=def-ncowan
#SBATCH --job-name=mm_g1c
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --output=/scratch/wasi/pcm/logs/g1c-%A_%a.out
#SBATCH --array=0-1
# G1 rerun: the first pass wrote no temperature because diagfi.def asked for
# "temp", which the generic physics does not write (it writes "temperature").
set -uo pipefail
TAGS=(g1inh g1noinh)
tag=${TAGS[$SLURM_ARRAY_TASK_ID]}
d=/scratch/wasi/pcm/runs/T265_Z100_${tag}_1d
source /scratch/wasi/pcm/trunk/LMDZ.COMMON/arch/arch-FIR.env
ulimit -s unlimited
export OMP_STACKSIZE=400M OMP_NUM_THREADS=1
cd "$d" || exit 1
echo "=== $tag ==="; cat diagfi.def
rm -f diagfi.nc
./rcm1d.e > lrun_g1 2>&1
echo "rc=$?"
grep -viE "^#[0-9]|^ *$" lrun_g1 | tail -8
ls -la *.nc
