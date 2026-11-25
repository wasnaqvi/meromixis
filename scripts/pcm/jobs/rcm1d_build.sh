#!/bin/bash
#SBATCH --account=def-ncowan
#SBATCH --job-name=mm_rcm1dbuild
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4G
#SBATCH --output=/scratch/%u/pcm/logs/rcm1dbuild-%j.out
#
# rcm1d with more levels:  sbatch --export=ALL,NLEV=160 rcm1d_build.sh
# At 40 levels the inhibited layer of a kappa ~ p column (3.9-6.5 bar in the emulator)
# is ~1.5 levels thick, too thin for the moist adjustment to act on.  160 levels give
# ~0.035 dex per level, ~6 across it.  The executable is named by its level count
# (rcm1d_<NLEV>_phygeneric_b20x25_seq.e), so the 40-level one is left in place.
set -uo pipefail
NLEV=${NLEV:-160}
P=/scratch/$USER/pcm
source $P/trunk/LMDZ.COMMON/arch/arch-FIR.env
module load perl/5.36.1
export PATH=$P/tools/modipsl_fcm/modipsl/modeles/LMDZ/tools/fcm/bin:$PATH
ulimit -s unlimited

cd $P/trunk/LMDZ.COMMON
./makelmdz_fcm -arch FIR -p generic -d $NLEV -b 20x25 rcm1d > $P/build_rcm1d_$NLEV.log 2>&1
rc=$?
echo "build rc=$rc"
grep -iE "Build command finished|Error [0-9]" $P/build_rcm1d_$NLEV.log | tail -3
ls -la bin/ | grep rcm1d
exit $rc
