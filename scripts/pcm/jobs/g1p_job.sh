#!/bin/bash
#SBATCH --account=def-ncowan
#SBATCH --job-name=mm_g1p
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --output=/scratch/%u/pcm/logs/g1p-%A_%a.out
#SBATCH --array=0-3
# G1 with grey kappa ~ p (scripts/pcm/make_greyp_1d.py): identical-start twins in a
# column that convects.  With constant kappa nothing convected and the on/off twins
# were bit-identical.  20 000 days (g1e ran 100 000 days in ~4.5 h).  The log also
# reports whether the corrk table was read and whether each adjustment ever fired.
# NLEV (sbatch --export=ALL,NLEV=160) picks the rcm1d build and the _L<NLEV> run
# directories; 40, the default, took ~1.5 h, and the time limit covers 160.
set -uo pipefail
TAGS=(g1pon_lidstart g1poff_lidstart g1pon_mixstart g1poff_mixstart)
NLEV=${NLEV:-40}
SUF=$([ "$NLEV" = 40 ] || echo "_L$NLEV")
tag=${TAGS[$SLURM_ARRAY_TASK_ID]}
d=/scratch/$USER/pcm/runs/T265_Z100_${tag}${SUF}_1d
source /scratch/$USER/pcm/trunk/LMDZ.COMMON/arch/arch-FIR.env
ulimit -s unlimited
export OMP_STACKSIZE=400M OMP_NUM_THREADS=1
cd "$d" || exit 1
echo "=== $tag, $NLEV levels ==="
ln -sfn /scratch/$USER/pcm/datagcm datagcm
ln -sf /scratch/$USER/pcm/trunk/LMDZ.COMMON/bin/rcm1d_${NLEV}_phygeneric_b20x25_seq.e rcm1d.e
grep -E "^(ichoice|ndt|moist_convection_inhibition|virtual_theta_correction)" rcm1d.def callphys.def
head -3 profile
rm -f diagfi.nc restartfi.nc profile.out
./rcm1d.e > lrun_g1 2>&1
echo "rc=$?"
grep -E "No initial profile|inifis: (moist_convection_inhibition|virtual_theta)|USING DEFAULTS : (ichoice|moist_convection|virtual_theta)" lrun_g1 | sort -u
grep -viE "^#[0-9]|^ *$" lrun_g1 | tail -8
grep -m3 -iE "corrk_gcm_IR|using no corrk|No absorption data" lrun_g1
ls -la *.nc
