#!/bin/bash
#SBATCH --account=def-ncowan
#SBATCH --job-name=mm_starts
#SBATCH --time=1:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4G
#SBATCH --output=/scratch/%u/pcm/logs/starts-%j.out
#
# Start files for the 3-D gate runs (G2, G3).  The lon-lat core cannot cold-start
# generic physics (iniacademic fills fields only for iflag_phys >= 2), so the first
# G2/G3 attempt stopped in dynetat0 on the missing start.nc.
#
#   1. build start2archive on the seed's grid (LMD reference_trappist1c, 64x48x50)
#   2. seed start.nc + startfi.nc -> start_archive.nc
#   3. per run: newstart (64x48x40, option 0) with the generated newstart.in, which
#      sets the planet constants, flattens the surface, puts p_s = p_bot, and loads
#      T and water from profile.in / profile_H2O_vap (the emulator column)
#
# Submit after newstart_build.sh has finished (two FCM builds must not share the tree).
set -uo pipefail
P=/scratch/$USER/pcm
SEED=$P/seed/reference_trappist1c_64x48x50_b32x36
RUNS="T265_Z100_g2dry T265_Z100_g3lid T265_Z100_g3noinh"

source $P/trunk/LMDZ.COMMON/arch/arch-FIR.env
module load perl/5.36.1
export PATH=$P/tools/modipsl_fcm/modipsl/modeles/LMDZ/tools/fcm/bin:$PATH
ulimit -s unlimited

echo "=============== BUILD start2archive (seed grid) ==============="
cd $P/trunk/LMDZ.COMMON
./makelmdz_fcm -arch FIR -p generic -d 64x48x50 -b 32x36 start2archive > $P/build_start2archive.log 2>&1
rc=$?; echo "build rc=$rc"
grep -iE "Build command finished|Error [0-9]" $P/build_start2archive.log | tail -3
S2A=$(ls $P/trunk/LMDZ.COMMON/bin/start2archive_64x48x50_phygeneric*seq.e 2>/dev/null | head -1)
NS=$(ls $P/trunk/LMDZ.COMMON/bin/newstart_64x48x40_phygeneric*seq.e 2>/dev/null | head -1)
echo "start2archive: ${S2A:-MISSING}"; echo "newstart: ${NS:-MISSING}"
[ -z "$S2A" ] || [ -z "$NS" ] && { echo "missing executable - stopping"; exit 1; }

echo "=============== start_archive.nc from the seed ==============="
cd $SEED || exit 1
ln -sfn $P/datagcm datadir
rm -f start_archive.nc
timeout 900 $S2A < /dev/null > lstart2archive 2>&1
echo "start2archive rc=$?"
tail -5 lstart2archive
[ -f start_archive.nc ] || { echo "no start_archive.nc - stopping"; exit 1; }

fail=0
for r in $RUNS; do
  echo "=============== newstart $r ==============="
  cd $P/runs/$r || { fail=1; continue; }
  rm -f start.nc startfi.nc restart.nc restartfi.nc
  cp $SEED/start_archive.nc .
  # timeout: at EOF some newstart prompts loop forever on a failed read
  timeout 600 $NS < newstart.in > lnewstart 2>&1
  echo "newstart rc=$?"
  grep -E "Could not find|problem reading|wrong input|not found|initialized using" lnewstart | head
  if [ -f restart.nc ] && [ -f restartfi.nc ]; then
    mv restart.nc start.nc; mv restartfi.nc startfi.nc
    rm -f start_archive.nc
    ls -la start.nc startfi.nc
  else
    echo "$r: newstart wrote no restart files"; tail -20 lnewstart; fail=1
  fi
done
exit $fail
