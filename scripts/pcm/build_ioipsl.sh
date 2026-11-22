#!/bin/bash
# Build IOIPSL on Fir from the LMD mirror.
# forge.ipsl.fr (the canonical SVN source used by install_ioipsl_*.bash) is
# unreachable from Fir, so we take the source snapshot LMD mirrors under
# ~lmdz/Distrib/IOIPSL/ and compile it directly with the arch-FIR flags.
set -euo pipefail
SRC=https://web.lmd.jussieu.fr/~lmdz/Distrib/IOIPSL/src_fairhead_v3c
ROOT=/scratch/$USER/pcm/trunk/IOIPSL
DL=/scratch/$USER/pcm/ioipsl_download          # untrusted download, kept separate

source /scratch/$USER/pcm/trunk/LMDZ.COMMON/arch/arch-FIR.env

rm -rf "$ROOT" "$DL"
mkdir -p "$DL" "$ROOT"/{lib,inc,bin,build}

cd "$DL"
for f in defprec.f90 stringop.f90 errioipsl.f90 mathelp.f90 calendar.f90 \
         getincom.f90 fliocom.f90 flincom.f90 histcom.f90 restcom.f90 ioipsl.f90; do
    curl -sS --max-time 90 -O "$SRC/$f"
done
ls -la "$DL" | tail -12

cd "$ROOT/build"
FFLAGS="-c -fdefault-real-8 -fdefault-double-8 -ffree-line-length-none -fno-align-commons -fallow-argument-mismatch -O2 -DNC_DOUBLE"
INC="$(nf-config --fflags)"
# dependency order: each module must be compiled after the modules it uses
for f in defprec stringop errioipsl mathelp calendar getincom fliocom flincom histcom restcom ioipsl; do
    echo "--- $f"
    mpif90 $FFLAGS $INC -J"$ROOT/inc" -I"$ROOT/inc" "$DL/$f.f90" -o "$f.o" || { echo "FAILED on $f"; exit 1; }
done
ar rcs "$ROOT/lib/libioipsl.a" *.o
echo "=== built ==="
ls -la "$ROOT/lib/" && ls "$ROOT/inc/" | head -20
