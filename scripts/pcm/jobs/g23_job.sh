#!/bin/bash
#SBATCH --account=def-ncowan
#SBATCH --job-name=mm_g23
#SBATCH --time=3:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=2G
#SBATCH --output=/scratch/$USER/pcm/logs/g23-%A_%a.out
#SBATCH --array=0-2
#
# Gates G2 (dry dynamics) and G3 (lid survival, with its no-inhibition control),
# double-grey, 64x48x40.  One array task per run.
set -uo pipefail
RUNS=(T265_Z100_g2dry T265_Z100_g3lid T265_Z100_g3noinh)
d=/scratch/$USER/pcm/runs/${RUNS[$SLURM_ARRAY_TASK_ID]}

source /scratch/$USER/pcm/trunk/LMDZ.COMMON/arch/arch-FIR.env
ulimit -s unlimited
export OMP_STACKSIZE=400M
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}

cd "$d" || exit 1
echo "=== $(basename $d) in $d ==="
ln -sf /scratch/$USER/pcm/datagcm datagcm
ln -sf /scratch/$USER/pcm/trunk/LMDZ.COMMON/bin/gcm_64x48x40_phygeneric_b20x25_para.e gcm.e
ls -la start.nc startfi.nc 2>/dev/null || echo "(no start files; relying on the model's own init)"

srun ./gcm.e > lrun_g23 2>&1
echo "gcm rc=$?"
grep -viE "^#[0-9]|^ *$" lrun_g23 | tail -25
echo "--- outputs ---"
ls -la *.nc 2>/dev/null
