#!/bin/bash
#SBATCH --account=def-ncowan
#SBATCH --job-name=mm_g1
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4G
#SBATCH --output=/scratch/%u/pcm/logs/g1-%j.out
#
# Gate G1: rebuild rcm1d (production, not -debug) and run the two column twins.
# Runs under Slurm so an SSH disconnect cannot kill it.
set -uo pipefail
mkdir -p /scratch/$USER/pcm/logs

source /scratch/$USER/pcm/trunk/LMDZ.COMMON/arch/arch-FIR.env
module load perl/5.36.1
export PATH=/scratch/$USER/pcm/tools/modipsl_fcm/modipsl/modeles/LMDZ/tools/fcm/bin:$PATH
ulimit -s unlimited
export OMP_STACKSIZE=400M OMP_NUM_THREADS=1

echo "=============== REBUILD rcm1d (no -debug) ==============="
cd /scratch/$USER/pcm/trunk/LMDZ.COMMON
./makelmdz_fcm -arch FIR -p generic -d 40 -b 20x25 rcm1d > /scratch/$USER/pcm/build_rcm1d.log 2>&1
rc=$?
echo "build rc=$rc"
grep -iE "Build command finished|Error [0-9]" /scratch/$USER/pcm/build_rcm1d.log | tail -3
ls -la bin/ | grep rcm1d
[ $rc -ne 0 ] && { echo "BUILD FAILED - stopping"; exit 1; }

echo
echo "=============== VERIFY the H2O tracer row ==============="
tail -2 /scratch/$USER/pcm/datagcm/table_tracers_condensable

for tag in g1inh g1noinh; do
  d=/scratch/$USER/pcm/runs/T265_Z100_${tag}_1d
  echo
  echo "=============== RUN $tag ==============="
  cd "$d" || { echo "missing $d"; continue; }
  ln -sf /scratch/$USER/pcm/trunk/LMDZ.COMMON/bin/rcm1d_40_phygeneric_b20x25_seq.e rcm1d.e
  rm -f diagfi.nc
  ./rcm1d.e > lrun_g1 2>&1
  echo "$tag rc=$?"
  grep -viE "^#[0-9]|^ *$" lrun_g1 | tail -12
  echo "--- outputs ---"
  ls -la *.nc 2>/dev/null
done
echo
echo "=============== DONE ==============="
