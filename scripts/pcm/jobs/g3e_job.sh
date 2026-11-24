#!/bin/bash
#SBATCH --account=def-ncowan
#SBATCH --job-name=mm_g3e
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G
#SBATCH --output=/scratch/%u/pcm/logs/g3e-%A_%a.out
#SBATCH --array=0-3
#
# G3, second pass: identical-start twins (scripts/pcm/make_twins.py), 500 days each
# (~5 h; the first pass ran 50 days in 25-33 min).  The start files are the ones
# make_starts.sh wrote for the old g3lid (lid column) and g3noinh (mixed column):
# newstart's input does not depend on the switch, so they are reused, after checking
# that the column they were made from is the one this run asks for.
set -uo pipefail
RUNS=(g3on_lidstart g3off_lidstart g3on_mixstart g3off_mixstart)
run=${RUNS[$SLURM_ARRAY_TASK_ID]}
case $run in
  *_lidstart) src=T265_Z100_g3lid ;;
  *_mixstart) src=T265_Z100_g3noinh ;;
esac
R=/scratch/$USER/pcm/runs
d=$R/T265_Z100_$run

source /scratch/$USER/pcm/trunk/LMDZ.COMMON/arch/arch-FIR.env
ulimit -s unlimited
export OMP_STACKSIZE=400M
# One thread per rank: see g23_job.sh (table_tracers_condensable on a fixed unit).
export OMP_NUM_THREADS=1

cd "$d" || exit 1
echo "=== $run in $d, start files from $src ==="
for f in profile.in profile_H2O_vap newstart.in; do
  cmp -s $f $R/$src/$f || { echo "$f differs from $src/$f: wrong start column"; exit 1; }
done
cp $R/$src/start.nc $R/$src/startfi.nc .
ln -sfn /scratch/$USER/pcm/datagcm datagcm
ln -sf /scratch/$USER/pcm/trunk/LMDZ.COMMON/bin/gcm_64x48x40_phygeneric_b20x25_para.e gcm.e
grep -E "^ *(nday|moist_convection_inhibition|virtual_theta_correction) *=" run.def callphys.def
rm -f lrun_g3 used_*.def diagfi.nc restart.nc restartfi.nc

srun ./gcm.e > lrun_g3 2>&1
echo "gcm rc=$?"
grep -iE "^ *(moist_convection_inhibition|virtual_theta_correction|nday) *=" used_callphys.def used_run.def
grep -viE "^#[0-9]|^ *$" lrun_g3 | tail -25
echo "--- outputs ---"
ls -la *.nc 2>/dev/null
