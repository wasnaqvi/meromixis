#!/bin/bash
#SBATCH --account=def-ncowan
#SBATCH --job-name=mm_g1d
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --output=/scratch/%u/pcm/logs/g1d-%A_%a.out
#SBATCH --array=0-1
# G1, third pass.  The previous pass started from rcm1d's isothermal column
# (ichoice=1, tref=397.5 K) with water only in the bottom layer.  Stable from the
# start, the column never convected, water never left the bottom layer, and the
# inhibition switch had nothing to act on: the twins differed by <0.1 K after
# 13 000 days.  The run directories now start from the emulator equilibrium
# (ichoice=8 with "profile" and "profile_H2O_vap"; deffiles.initial_column).
set -uo pipefail
TAGS=(g1inh g1noinh)
tag=${TAGS[$SLURM_ARRAY_TASK_ID]}
d=/scratch/$USER/pcm/runs/T265_Z100_${tag}_1d
source /scratch/$USER/pcm/trunk/LMDZ.COMMON/arch/arch-FIR.env
ulimit -s unlimited
export OMP_STACKSIZE=400M OMP_NUM_THREADS=1
cd "$d" || exit 1
echo "=== $tag ==="
grep -E "^(ichoice|moist_convection_inhibition|enertest)" rcm1d.def callphys.def
rm -f diagfi.nc restartfi.nc profile.out
./rcm1d.e > lrun_g1 2>&1
echo "rc=$?"
grep -E "No initial profile|USING DEFAULTS : (ichoice|moist_convection|virtual_theta)" lrun_g1
grep -viE "^#[0-9]|^ *$" lrun_g1 | tail -8
ls -la *.nc
