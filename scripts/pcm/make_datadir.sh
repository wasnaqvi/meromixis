#!/bin/bash
# Build a minimal datagcm for double-grey Generic PCM runs on Fir.
#
# No public corrk table has variable water above x_H2O = 1e-2 in an H2/He
# background (the LMD repository's LMDZ_Jupiter_water_var stops there), and we need
# 0.089 at 100x solar and 0.30 at 300x.  So the correlated-k spine is not runnable.
#
# What IS runnable is the model's double-grey mode, which is exactly the `grey` gate
# run in experiments/steam_rising.yaml.  rad_correlatedk_read_opacity_tables.F90
# zeroes the gaseous k-coefficients when corrkdir starts with "null" (line ~409) and
# then applies the constant kappa_IR / kappa_VI when graybody = .true.  It still
# reads the axis and band-structure files unconditionally, so we assemble a corrk
# directory that has those and no k-data, borrowing the band structure of the
# Jupiter table so the bands are physically sensible.
set -euo pipefail
B=https://web.lmd.jussieu.fr/~lmdz/planets/generic/datagcm/corrk_data/LMDZ_Jupiter_water_var
D=/scratch/$USER/pcm/datagcm/corrk_data/null_H2He_20x25

mkdir -p "$D/20x25"
cd "$D"
for f in T.dat p.dat Q.dat g.dat; do curl -sS --max-time 90 -O "$B/$f"; done
cd "$D/20x25"
for f in narrowbands_IR.in narrowbands_VI.in; do curl -sS --max-time 90 -O "$B/20x25/$f"; done

echo "=== assembled ==="
find /scratch/$USER/pcm/datagcm -type f | sed "s|/scratch/$USER/pcm/datagcm/||"
echo
echo "--- band counts (must match the -b flag used to build) ---"
echo "IR bands: $(head -1 "$D/20x25/narrowbands_IR.in")"
echo "VI bands: $(head -1 "$D/20x25/narrowbands_VI.in")"
echo "g points: $(head -1 "$D/g.dat")"
echo "T nodes : $(head -1 "$D/T.dat")"
echo "p nodes : $(head -1 "$D/p.dat")"
