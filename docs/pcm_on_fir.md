# Running the sweep on DRAC Fir

A short runbook. The project primer has the background; this lists the commands in order.

## Once

```bash
cd /project/def-ncowan/$USER
svn checkout https://svn.lmd.jussieu.fr/Planeto/trunk --depth empty && cd trunk
svn update LMDZ.COMMON LMDZ.GENERIC
svn info | grep Revision            # record it; this repo's templates were checked at r4146
```

Arch files (`LMDZ.COMMON/arch/arch-FIR.{env,fcm,path}`) are in the primer. Build IOIPSL,
then the benchmark, then:

```bash
cd LMDZ.COMMON
./makelmdz_fcm -arch FIR -p generic -d 40 -b 38x36 rcm1d
./makelmdz_fcm -arch FIR -p generic -d 64x48x40 -b 38x36 -parallel mpi_omp gcm
```

The physics directory was renamed `phygeneric` in the trunk; if `-p generic` is refused on
your revision use `-p std`. One executable per grid, band count and tracer-aerosol count
(`-t 1` for the `clouds` variants).

## Per campaign

```bash
pip install -e /path/to/meromixis
meromixis pcm-init experiments/steam_rising.yaml /scratch/$USER/steam-rising
column -t -s $'\t' /scratch/$USER/steam-rising/runs.tsv | less -S
```

Each run directory holds `callphys.def`, `run.def` (or `rcm1d.def`), `gases.def`,
`traceur.def`, `z2sig.def`, `diagfi.def`, `job.slurm`, `planet_constants.txt` and
`manifest.json`. Before the first submission:

1. Diff `callphys.def` and `run.def` against the nearest `LMDZ.GENERIC/deftank` example
   at your revision (`docs/VERIFY.md`, items P1-P7).
2. Set `corrkdir` in the YAML to a real H2-He table with variable water and regenerate.
3. Link the executable as `gcm.e` (or `rcm1d.e`) and `datadir` into each directory.
4. Make start files: run the `_1d` twin, then `newstart` with the values in
   `planet_constants.txt`, initialising T(p) and q(p) from the column. Runs with a
   `branch_from` entry copy `segments/restart<N>.nc` and `restartfi<N>.nc` of the parent to
   `start.nc` and `startfi.nc` instead, and change nothing but `Fat1AU`.

## Running

```bash
cd /scratch/$USER/steam-rising/T300_Z100_fid && sbatch job.slurm
tail -n 5 spinup.log                # one line per finished segment
```

The job runs one segment, checks for `Everything is cool`, archives `diagfi` and restarts
under `segments/`, promotes the restart files and resubmits itself until `num_run` reaches
`NUMMAX`. A failed segment stops the chain and leaves the last good restart in place.
Copy `segments/` to `/project` periodically; `/scratch` is purged.

## Analysing

```bash
meromixis diag segments/diagfi{37,38,39,40}.nc --gravity 12.41 --metallicity 100 \
    --varmap q_vap=H2O_vap q_ice=H2O_ice --out summary_T300_Z100.nc
```

`--varmap` renames tracers for the generic route; check names once with `ncdump -h`.
If `w` is written in Pa s⁻¹ add `--w-is-omega`.
