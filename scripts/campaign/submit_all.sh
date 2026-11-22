#!/bin/bash
# Submit the gate-0 scoping campaign on Fir.
#   bash scripts/campaign/submit_all.sh
# Jobs are independent; C is tiny and runs alone, A and B are array jobs.
set -euo pipefail

ACCOUNT=def-ncowan
RUN=${MEROMIXIS_RUN:-/scratch/wasi/meromixis/gate0}
REPO=${MEROMIXIS_REPO:-$HOME/meromixis}
mkdir -p "$RUN"/{job_a,job_b,job_c,logs}

common=(--account=$ACCOUNT --nodes=1 --ntasks=1 --cpus-per-task=1 --mem-per-cpu=2G)

jid_a=$(sbatch --parsable "${common[@]}" --job-name=mm_jobA --time=0:30:00 \
  --array=0-63 --output="$RUN/logs/jobA-%A_%a.out" \
  --wrap="$REPO/scripts/campaign/run_task.sh \
          $REPO/scripts/campaign/job_a_transition.py \$SLURM_ARRAY_TASK_ID 64 $RUN/job_a")
echo "job A (transition ensemble, 512 draws / 64 tasks): $jid_a"

jid_b=$(sbatch --parsable "${common[@]}" --job-name=mm_jobB --time=0:30:00 \
  --array=0-15 --output="$RUN/logs/jobB-%A_%a.out" \
  --wrap="$REPO/scripts/campaign/run_task.sh \
          $REPO/scripts/campaign/job_b_barrier.py \$SLURM_ARRAY_TASK_ID 16 $RUN/job_b")
echo "job B (barrier kill test, 256 draws / 16 tasks): $jid_b"

jid_c=$(sbatch --parsable "${common[@]}" --job-name=mm_jobC --time=1:00:00 \
  --output="$RUN/logs/jobC-%j.out" \
  --wrap="$REPO/scripts/campaign/run_task.sh \
          $REPO/scripts/campaign/job_c_checks.py $RUN/job_c")
echo "job C (VERIFY checks): $jid_c"

echo "$jid_a $jid_b $jid_c" > "$RUN/jobids.txt"
echo "run dir: $RUN"
