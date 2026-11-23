#!/bin/bash
#SBATCH --account=def-ncowan
#SBATCH --job-name=mm_nsbuild
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4G
#SBATCH --output=/scratch/%u/pcm/logs/nsbuild-%j.out
#
# Build newstart on the 3-D grid.  The lon-lat core cannot cold-start generic
# physics: with read_start=.false. iniacademic only fills fields when
# iflag_phys >= 2, so every 3-D run needs start.nc/startfi.nc made by newstart.
set -uo pipefail
mkdir -p /scratch/$USER/pcm/logs

source /scratch/$USER/pcm/trunk/LMDZ.COMMON/arch/arch-FIR.env
module load perl/5.36.1
export PATH=/scratch/$USER/pcm/tools/modipsl_fcm/modipsl/modeles/LMDZ/tools/fcm/bin:$PATH
ulimit -s unlimited

cd /scratch/$USER/pcm/trunk/LMDZ.COMMON
./makelmdz_fcm -arch FIR -p generic -d 64x48x40 -b 20x25 newstart > /scratch/$USER/pcm/build_newstart.log 2>&1
rc=$?
echo "build rc=$rc"
grep -iE "Build command finished|Error [0-9]" /scratch/$USER/pcm/build_newstart.log | tail -3
ls -la bin/ | grep newstart
exit $rc
