#!/bin/bash
#SBATCH --account=def-ncowan
#SBATCH --job-name=mm_g1e
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --output=/scratch/%u/pcm/logs/g1e-%A_%a.out
#SBATCH --array=0-3
# G1, fourth pass: identical-start twins (scripts/pcm/make_twins.py).  In g1d each twin
# started from the emulator column for its own switch, and the 8 K bottom difference
# it reported was there at day 50 and barely moved by day 19 200.  Here each on/off
# pair shares one start column, from both columns, for 100 000 days (~6.5 h; g1d ran
# 19 200 days in 73 min).
set -uo pipefail
TAGS=(g1on_lidstart g1off_lidstart g1on_mixstart g1off_mixstart)
tag=${TAGS[$SLURM_ARRAY_TASK_ID]}
d=/scratch/$USER/pcm/runs/T265_Z100_${tag}_1d
source /scratch/$USER/pcm/trunk/LMDZ.COMMON/arch/arch-FIR.env
ulimit -s unlimited
export OMP_STACKSIZE=400M OMP_NUM_THREADS=1
cd "$d" || exit 1
echo "=== $tag ==="
ln -sfn /scratch/$USER/pcm/datagcm datagcm
ln -sf /scratch/$USER/pcm/trunk/LMDZ.COMMON/bin/rcm1d_40_phygeneric_b20x25_seq.e rcm1d.e
grep -E "^(ichoice|ndt|moist_convection_inhibition|virtual_theta_correction)" rcm1d.def callphys.def
head -3 profile
rm -f diagfi.nc restartfi.nc profile.out
./rcm1d.e > lrun_g1 2>&1
echo "rc=$?"
grep -E "No initial profile|inifis: (moist_convection_inhibition|virtual_theta)|USING DEFAULTS : (ichoice|moist_convection|virtual_theta)" lrun_g1 | sort -u
grep -viE "^#[0-9]|^ *$" lrun_g1 | tail -8
ls -la *.nc
