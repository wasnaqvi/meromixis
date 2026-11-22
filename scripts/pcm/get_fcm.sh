#!/bin/bash
# Fetch FCM 1.x for LMDZ on Fir.
#
# LMDZ needs FCM 1.x.  The current metomi/fcm on GitHub is FCM 2.x and rejects
# LMDZ's bld.cfg ("bld::tool::fflags::muphy: invalid sub-package in declaration",
# "dir::bin: label not recognised").  FCM 1.x is not tagged there and the Planeto
# trunk has no tools/ directory, but FCM 1.x ships inside the IPSL modipsl bundle,
# so stream the bundle and extract only that subtree.
#
# Then fix the shebangs: Fir's /usr/bin/perl 5.32 is stripped and has no FindBin,
# while fcm and fcm_internal both hardcode #!/usr/bin/perl.  Wrapping just `fcm`
# is not enough because the generated Makefile invokes fcm_internal directly.
set -euo pipefail
DEST=${1:-/scratch/$USER/pcm/tools/modipsl_fcm}
URL=https://lmdz.lmd.jussieu.fr/pub/src_archives/unstable/modipsl.20251128.trunk.tar.gz

mkdir -p "$DEST" && cd "$DEST"
curl -sS --max-time 900 "$URL" | tar -xz --wildcards "*/tools/fcm/*"

FCM=$(find "$DEST" -type d -name fcm -path "*/tools/*" | head -1)
[ -n "$FCM" ] || { echo "fcm not found in bundle" >&2; exit 1; }
sed -i "1s|^#!/usr/bin/perl|#!/usr/bin/env perl|" "$FCM"/bin/*
echo "FCM 1.x at: $FCM/bin"
echo "Use with:  module load perl/5.36.1 && export PATH=$FCM/bin:\$PATH"
