# STAR-Fusion WDL Pipeline — User Guide

This guide describes how to detect gene fusions from paired-end RNA-seq data using
the STAR-Fusion WDL pipeline (repository:
[`AmyOlex/STAR-Fusion-WDL-hg38v22`](https://github.com/AmyOlex/STAR-Fusion-WDL-hg38v22)).
It consolidates the WDL pipeline reference, the HPC execution walkthrough, and the
batch execution SOP into a single end-to-end reference.

> ⚠️ **About the paths in this guide:** The example commands below are filled in with the
> real paths from a run that completed successfully on Apollo — the
> `COAD_PAAD_FusionAnalysis4Yang_260504` project (see
> [Section 1.1](#11-a-working-example-on-apollo)) — so you have a concrete, known-good
> reference to follow. **When you set up a new project you must update these paths to your
> own.** The ones that change per project/sample are your run/config directory (created
> under `.../config_wdl/`), your batch `output_dir` (under `.../WDL_STARFusion_Results/`),
> the `sample_id`, and the input FASTQ (or junction) file paths. The reference genome and
> the pipeline repository paths normally stay the same across projects unless you relocate
> them (see the notes in Sections 5 and 8.3).

---

## 1. Introduction

The pipeline wraps **STAR-Fusion v1.10** (with optional FusionInspector validation and
coding-effect annotation) in WDL so that the same workflow logic can run in two
environments:

- **Terra / Google Cloud** — using the cloud wrapper, with the CTAT genome library
  supplied as a tar.gz and extracted at runtime.
- **Local HPC (SLURM + Cromwell + Singularity)** — using the HPC wrapper, with a
  pre-extracted CTAT genome library on the local filesystem.

**A note on the terminology, if this is new to you:** **WDL** (Workflow Description
Language) is a plain-text format for describing an analysis as a series of steps — what
commands to run, what inputs each step needs, and how the steps connect. Writing the
pipeline in WDL means the same recipe can run unchanged in different places. **Cromwell**
is the program that actually reads a WDL file and executes it: it hands each step to the
compute system, tracks which steps have finished, and gathers the results. On the HPC,
Cromwell submits jobs to **SLURM** (the cluster's job scheduler) and runs each step inside
a **Singularity** container (a self-contained software environment that bundles STAR-Fusion
and its dependencies, so you don't install them yourself). On Terra, this same role is
handled for you by the platform. In short: you describe *what* to run in a WDL and its
inputs file, and Cromwell takes care of *running* it.

On the local HPC you can run the pipeline in two modes:

- **Single-sample** — one FASTQ pair (or one junction file) per run.
- **Batch** — many samples in parallel, driven by a CSV sample sheet. Batch mode uses a
  scatter workflow and a few extra preparation steps, described in
  [Section 8](#8-running-on-the-hpc--batch).

The two wrappers share the same underlying STAR-Fusion task, so results are equivalent
across environments; only the genome input, file paths, and job-submission mechanics
differ.

### 1.1 A working example on Apollo

A complete run that finished successfully on Apollo is kept as a reference here:

```
/lustre/home/harrell_lab/bulkRNASeq/config_wdl/COAD_PAAD_FusionAnalysis4Yang_260504
```

If anything in the setup below is unclear, look at this directory for a concrete, known-good
example of the `cromwell.conf`, inputs/batch files, and `run_command.sh` for a real run.

**Directory conventions on Apollo:**

- **Config and run directories** live under
  `/lustre/home/harrell_lab/bulkRNASeq/config_wdl/` — create a new subdirectory here for
  each analysis (as in the example above), so all WDL configuration stays in one place.
- **Results** are saved under
  `/lustre/home/harrell_lab/bulkRNASeq/WDL_STARFusion_Results/` — point your batch
  `output_dir` (Section 8.2) into this directory so gathered results collect in a
  consistent location.

---

## 2. Scope and Important Caveats

**This pipeline calls fusions using a single tool: STAR-Fusion.** It does not run
multiple fusion callers. Many published fusion analyses combine two or more callers
(e.g., STAR-Fusion, Arriba, FusionCatcher) and take a consensus. If your project
requires multi-caller consensus, that is outside the scope of this pipeline and must be
added separately.

**The output is a set of raw fusion predictions, not a finished result set.** STAR-Fusion
predictions still require downstream filtering, annotation, and analysis before they can
be interpreted — for example, applying read-support and FFPM thresholds, removing likely
artifacts and same-gene-family events, prioritizing by known-fusion databases, and
performing any cohort-level or expression analysis. See
[Section 13](#13-downstream-analysis-and-filtering) for where that downstream work lives
and important guidance on reusing it.

---

## 3. The WDL Pipeline

### 3.1 WDL files

The repository contains three WDL files:

| File | Purpose |
|------|---------|
| `star_fusion_workflow.wdl` | Sub-workflow containing the STAR-Fusion task logic. Shared by both wrappers. |
| `star_fusion_hg38_wf.wdl` | **Terra / Cloud wrapper.** Accepts `genome_plug_n_play_tar_gz` (tar.gz, default points to a GCS-hosted build). Extracts the genome at runtime. |
| `star_fusion_hg38_hpc_wf.wdl` | **HPC wrapper.** Accepts `local_genome_dir` (a String path to a pre-extracted CTAT genome library). No runtime extraction — saves disk space and time. |

The sub-workflow supports both input modes. The two wrappers provide clean entry points
for each environment without user-level flags or conditional logic. **On HPC, always use
`star_fusion_hg38_hpc_wf.wdl`.**

- **Docker image:** `trinityctat/starfusion:latest`
- **Default resources:** 12 CPUs, 100 GB memory (STAR-Fusion typically needs 50–100 GB)
- **CTAT reference library:** `GRCh38_gencode_v22_CTAT_lib_Mar012021` (plug-n-play)

### 3.2 Input modes

**Mode 1 — FASTQ input (standard).** Provide either:

- `left_fq` + `right_fq` — paired FASTQ files, or
- `fastq_pair_tar_gz` — a tarball of paired FASTQs

**Mode 2 — Junction file re-run.** Provide `input_chimeric_junction`
(`Chimeric.out.junction` or `.gz`) from a prior STAR-Fusion run. STAR-Fusion is invoked
with `-J`, skipping STAR alignment entirely. This is useful for re-running FusionInspector
validation or coding-effect annotation without repeating the expensive alignment step.
The `Chimeric.out.junction.gz` produced by this pipeline can be fed directly back in.

**Genome input depends on environment:**

- **Terra/Cloud:** `genome_plug_n_play_tar_gz` (extracted at runtime).
- **HPC:** `local_genome_dir` — full path to a pre-extracted `ctat_genome_lib_build_dir`.

---

## 4. Inputs and Outputs Reference

### 4.1 HPC wrapper inputs (`star_fusion_hg38_hpc_wf.wdl`)

| Parameter | Type | Required | Default | Description |
|---|---|---|---|---|
| `sample_id` | String | Yes | | Sample identifier |
| `local_genome_dir` | String | Yes | | Path to pre-extracted `ctat_genome_lib_build_dir` |
| `left_fq` | File? | No* | | R1 FASTQ file |
| `right_fq` | File? | No* | | R2 FASTQ file |
| `fastq_pair_tar_gz` | File? | No* | | Tarball of paired FASTQs |
| `input_chimeric_junction` | File? | No* | | `Chimeric.out.junction(.gz)` from a prior run |
| `fusion_inspector` | String? | No | | `inspect` or `validate` |
| `examine_coding_effect` | Boolean | No | `false` | Annotate coding effects |
| `coord_sort_bam` | Boolean | No | `false` | Coordinate-sort the output BAM |
| `min_FFPM` | Float | No | `0.1` | Minimum fusion fragments per million |
| `docker` | String | No | `trinityctat/starfusion:latest` | Docker image |
| `num_cpu` | Int | No | `12` | Number of CPUs |
| `memory` | String | No | `50G` | Memory allocation |

\* One of `left_fq` + `right_fq`, `fastq_pair_tar_gz`, or `input_chimeric_junction` must
be provided.

The Terra/Cloud wrapper (`star_fusion_hg38_wf.wdl`) takes the same parameters, except it
uses `genome_plug_n_play_tar_gz` (File) instead of `local_genome_dir`, and adds
cloud-specific parameters: `preemptible` (default `2`), `use_ssd` (default `true`), and
`extra_disk_space` (default `10` GB). These cloud parameters are not present in the HPC
wrapper.

### 4.2 Outputs

| Output | Type | Description |
|---|---|---|
| `fusion_predictions` | File | STAR-Fusion predictions (gzipped TSV) |
| `fusion_predictions_abridged` | File | Abridged predictions (gzipped TSV) |
| `junction` | File? | `Chimeric.out.junction` (gzipped) — reusable for re-runs |
| `bam` | File? | Aligned BAM |
| `bai` | File? | BAM index (if coordinate-sorted) |
| `sj` | File? | `SJ.out.tab` splice-junction file (gzipped) |
| `star_log_final` | File? | STAR alignment log |
| `coding_effect` | File? | Coding-effect annotations (if `examine_coding_effect` = true) |
| `extract_fusion_reads` | Array[File]? | Fusion evidence reads |
| `fusion_inspector_validate_fusions_abridged` | File? | FusionInspector validate results |
| `fusion_inspector_validate_web` | File? | FusionInspector validate web report |
| `fusion_inspector_inspect_fusions_abridged` | File? | FusionInspector inspect results |
| `fusion_inspector_inspect_web` | File? | FusionInspector inspect web report |

Optional outputs (`File?`) are produced only when the relevant run mode is used (e.g., BAM
and junction files are produced only during full alignment runs, not junction re-runs).

---

## 5. Prerequisites (HPC)

Before setting up a run on the cluster, ensure the following are in place:

- **Java 17+** available on your `PATH`.
- **Cromwell 91** installed (e.g., at `~/.cromwell/lib/cromwell-91.jar`) with a working
  SLURM configuration. See the Cromwell HPC Installation Guide for setup and the
  hello-world test.
- **Singularity or Apptainer** available via `module load` on compute nodes.
- **CTAT genome library** downloaded and pre-extracted to a shared location.

Download and extract the CTAT library once (note: `gsutil` is not available on Apollo, so
use `wget`):

```bash
mkdir -p /lustre/home/harrell_lab/references/ctat   # current shared location on Apollo
cd /lustre/home/harrell_lab/references/ctat

wget https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/__genome_libs_StarFv1.10/GRCh38_gencode_v22_CTAT_lib_Mar012021.plug-n-play.tar.gz
tar xzf GRCh38_gencode_v22_CTAT_lib_Mar012021.plug-n-play.tar.gz
```

Verify the directory you will pass as `local_genome_dir`:

```bash
ls .../GRCh38_gencode_v22_CTAT_lib_Mar012021.plug-n-play/ctat_genome_lib_build_dir/
```

Once extracted and verified, you can delete the tar.gz to reclaim ~30 GB.

> **Where the reference lives (and using it in another project):** On the cluster the CTAT
> library currently resides under
> `/lustre/home/harrell_lab/references/ctat/GRCh38_gencode_v22_CTAT_lib_Mar012021.plug-n-play/ctat_genome_lib_build_dir`,
> and that is the path the `local_genome_dir` value in the inputs files below points to. The
> reference is generic and can be shared across projects, so it does not need to be
> re-downloaded per run. **If you are running this pipeline for a different project or lab,
> consider placing the extracted library in a more neutral, shared location** (rather than
> under a specific lab's directory) and update `local_genome_dir` in your inputs files
> accordingly.

**Clone the pipeline and switch to the HPC branch:**

```bash
cd /lustre/home/harrell_lab/src            # example path
git clone https://github.com/AmyOlex/STAR-Fusion-WDL-hg38v22.git
cd STAR-Fusion-WDL-hg38v22
git checkout hpc-local-genome

ls *.wdl
#   star_fusion_hg38_hpc_wf.wdl   (HPC wrapper — use this one)
#   star_fusion_hg38_wf.wdl       (Terra/Cloud wrapper)
#   star_fusion_workflow.wdl      (shared sub-workflow)
```

> **Note on FASTQ symlinks:** Symlinks are fine as long as the *real* files (check with
> `readlink -f <symlink>`) also resolve to a path under the bind-mounted directory.

---

## 6. Run Setup (HPC)

Create a fresh run directory for each analysis. The Cromwell execution tree, logs, and
config all live here.

```bash
mkdir -p /lustre/home/harrell_lab/bulkRNASeq/config_wdl/COAD_PAAD_FusionAnalysis4Yang_260504
cd /lustre/home/harrell_lab/bulkRNASeq/config_wdl/COAD_PAAD_FusionAnalysis4Yang_260504
```

### 6.1 The Cromwell configuration (`cromwell.conf`)

Create `cromwell.conf` in your run directory. The following SLURM + Singularity backend
works on Apollo and can be adapted for other clusters:

> ⚠️ **Before you copy this:** the `submit-docker` block below contains a hardcoded bind
> path, `--bind /lustre/home/harrell_lab:/lustre/home/harrell_lab`, which is an **example
> specific to the Apollo lab directory.** You must replace both halves of it with the base
> path where *your* data, references, and run directories live (see the "Customize these
> for your cluster" table just below the block). If this path doesn't cover your files, the
> container won't be able to see them and the run will fail.

```hocon
include required(classpath("application"))

backend {
  default = "SLURM"
  providers {
    SLURM {
      actor-factory = "cromwell.backend.impl.sfs.config.ConfigBackendLifecycleActorFactory"
      config {
        runtime-attributes = """
          Int cpu = 12
          String memory = "100G"
          String? docker
          String? slurm_partition = "cpu-small"
          String? time_limit = "24:00:00"
        """

        submit = """
          sbatch \
            --job-name=${job_name} \
            --partition=${slurm_partition} \
            --cpus-per-task=${cpu} \
            --mem=${memory} \
            --time=${time_limit} \
            --output=${cwd}/execution/main_stdout \
            --error=${cwd}/execution/main_stderr \
            /bin/bash ${script}
        """

        submit-docker = """
          echo '#!/bin/bash
          module load singularity
          singularity exec --bind ${cwd}:${docker_cwd} --bind /lustre/home/harrell_lab:/lustre/home/harrell_lab docker://${docker} ${job_shell} ${docker_script}' > ${cwd}/execution/submit_script.sh
          chmod 755 ${cwd}/execution/submit_script.sh
          sbatch \
            --job-name=${job_name} \
            --partition=${slurm_partition} \
            --cpus-per-task=${cpu} \
            --mem=${memory} \
            --time=${time_limit} \
            --output=${cwd}/execution/docker_stdout \
            --error=${cwd}/execution/docker_stderr \
            ${cwd}/execution/submit_script.sh
        """

        kill = "scancel ${job_id}"
        check-alive = "squeue -j ${job_id}"
        job-id-regex = "(\\d+)"

        root = "cromwell-executions"

        # Use symlinks instead of copying input files to save disk space
        filesystems {
          local {
            localization: [ "soft-link", "copy" ]
          }
        }
      }
    }
  }
}
```

Customize these for your cluster:

| What to change | How to find it |
|---|---|
| Partition (`slurm_partition`) | `sinfo -s`, or `sinfo -o "%P %l %c %m %a"` to find high-memory partitions |
| Singularity/Apptainer module name (`module load singularity`) | `module avail singularity` or `module avail apptainer` |
| Bind path (`--bind /base/path:/base/path`) | Must cover every path referenced in your inputs JSON (data, references, run dir) |
| `memory` / `cpu` / `time_limit` | Adjust to your partition limits and expected runtime (a full run can take several hours) |

> **File localization:** The `filesystems` block tells Cromwell to symlink input files
> instead of copying them, falling back to copy only if symlinks fail. This is essential
> on storage-constrained clusters.

> **Bind paths:** `--bind /base/path:/base/path` makes everything under that directory
> visible inside the container at the same path, so input paths in `inputs.json` resolve
> without translation. Each additional location needs its own `--bind` entry (e.g.,
> `--bind /path1:/path1 --bind /path2:/path2`). Symlinked FASTQs work only if the real
> files also live under a bound path.

### 6.2 A `run_command.sh` script

Record the exact execution command so the run is documented and easy to repeat:

```bash
cat > run_command.sh << 'EOF'
java -Dconfig.file=/lustre/home/harrell_lab/bulkRNASeq/config_wdl/COAD_PAAD_FusionAnalysis4Yang_260504/cromwell.conf \
  -jar ~/.cromwell/lib/cromwell-91.jar \
  run /lustre/home/harrell_lab/src/STAR-Fusion-WDL-hg38v22/star_fusion_hg38_hpc_wf.wdl \
  --inputs inputs.json 2>&1 | tee cromwell_run.log
EOF
```

For batch runs, point this at the batch WDL and `batch_inputs.json` instead (see
[Section 8](#8-running-on-the-hpc--batch)).

---

## 7. Running on the HPC — Single Sample

### 7.1 Create `inputs.json`

**Full FASTQ-to-fusion run:**

```json
{
  "star_fusion_hg38_wf.sample_id": "VCU-PC-124_9017",
  "star_fusion_hg38_wf.left_fq": "/path/to/sample_R1.fq.gz",
  "star_fusion_hg38_wf.right_fq": "/path/to/sample_R2.fq.gz",
  "star_fusion_hg38_wf.local_genome_dir": "/lustre/home/harrell_lab/references/ctat/GRCh38_gencode_v22_CTAT_lib_Mar012021.plug-n-play/ctat_genome_lib_build_dir",
  "star_fusion_hg38_wf.fusion_inspector": "validate",
  "star_fusion_hg38_wf.examine_coding_effect": true,
  "star_fusion_hg38_wf.memory": "100G",
  "star_fusion_hg38_wf.docker": "trinityctat/starfusion:latest"
}
```

**Junction re-run (skip alignment):** replace the FASTQ entries with a junction file.

```json
{
  "star_fusion_hg38_wf.sample_id": "VCU-PC-124_9017",
  "star_fusion_hg38_wf.input_chimeric_junction": "/path/to/Chimeric.out.junction",
  "star_fusion_hg38_wf.local_genome_dir": "/lustre/home/harrell_lab/references/ctat/GRCh38_gencode_v22_CTAT_lib_Mar012021.plug-n-play/ctat_genome_lib_build_dir",
  "star_fusion_hg38_wf.fusion_inspector": "validate",
  "star_fusion_hg38_wf.examine_coding_effect": true,
  "star_fusion_hg38_wf.memory": "100G",
  "star_fusion_hg38_wf.docker": "trinityctat/starfusion:latest"
}
```

All paths must be absolute and accessible from compute nodes via the `cromwell.conf` bind
mount.

### 7.2 Run

Cromwell runs as a foreground Java process that must stay alive for the whole pipeline.
If your SSH session drops, Cromwell dies and won't collect results. **Always run inside a
`tmux` session.**

```bash
tmux new -s starfusion
cd /lustre/home/harrell_lab/bulkRNASeq/config_wdl/COAD_PAAD_FusionAnalysis4Yang_260504
source run_command.sh

# Detach: Ctrl+B, then D      Reattach: tmux attach -t starfusion
```

| Action | Command |
|--------|---------|
| New session | `tmux new -s <name>` |
| Detach | `Ctrl+B`, then `D` |
| Reattach | `tmux attach -t <name>` |
| List sessions | `tmux ls` |
| Kill session | `tmux kill-session -t <name>` |

---

## 8. Running on the HPC — Batch

Batch mode runs multiple samples in parallel using the WDL `scatter` block. Cromwell
submits each sample as a separate SLURM job, and a gather task collects results into a
clean output directory. Use the batch WDL (`star_fusion_hg38_batch_hpc_wf.wdl`).

### 8.1 Prepare a CSV sample sheet

One row per sample. Required columns: `sample_id`, `left_fq`, `right_fq`.

```csv
sample_id,left_fq,right_fq
VCU-PC-124_9017,/lustre/home/harrell_lab/bulkRNASeq/raw/VCU-PC-124_9017_R1.fq.gz,/lustre/home/harrell_lab/bulkRNASeq/raw/VCU-PC-124_9017_R2.fq.gz
VCU-PC-125_9018,/lustre/home/harrell_lab/bulkRNASeq/raw/VCU-PC-125_9018_R1.fq.gz,/lustre/home/harrell_lab/bulkRNASeq/raw/VCU-PC-125_9018_R2.fq.gz
```

### 8.2 Edit `batch_config.json`

Parameters shared across all samples:

```json
{
  "star_fusion_hg38_batch_hpc_wf.local_genome_dir": "/lustre/home/harrell_lab/references/ctat/GRCh38_gencode_v22_CTAT_lib_Mar012021.plug-n-play/ctat_genome_lib_build_dir",
  "star_fusion_hg38_batch_hpc_wf.fusion_inspector": "validate",
  "star_fusion_hg38_batch_hpc_wf.examine_coding_effect": true,
  "star_fusion_hg38_batch_hpc_wf.memory": "100G",
  "star_fusion_hg38_batch_hpc_wf.docker": "trinityctat/starfusion:latest",
  "star_fusion_hg38_batch_hpc_wf.output_dir": "/lustre/home/harrell_lab/bulkRNASeq/WDL_STARFusion_Results/COAD_PAAD_FusionAnalysis4Yang_260504",
  "star_fusion_hg38_batch_hpc_wf.keep_bam_and_junction": false
}
```

- `output_dir` is where gathered results go. The gather task creates one subdirectory per
  sample:

  ```
  /lustre/home/harrell_lab/bulkRNASeq/WDL_STARFusion_Results/COAD_PAAD_FusionAnalysis4Yang_260504/
    VCU-PC-124_9017/
      VCU-PC-124_9017.star-fusion.fusion_predictions.tsv.gz
      VCU-PC-124_9017.star-fusion.fusion_predictions.abridged.tsv.gz
      ...
    VCU-PC-125_9018/
      ...
  ```

- Set `keep_bam_and_junction` to `true` to include BAM, junction, and SJ files in the
  gathered output. These are large and excluded by default.

### 8.3 Generate the batch inputs JSON

Use the helper script to convert the sample sheet + shared config into a combined inputs
file:

```bash
python /lustre/home/harrell_lab/src/STAR-Fusion-WDL-hg38v22/example_hpc_run/csv_to_batch_inputs.py \
  sample_sheet.csv \
  --config batch_config.json \
  --output batch_inputs.json
```

> **Where the scripts live (and using them in another project):** This helper script ships
> inside the cloned pipeline repository, which currently resides at
> `/lustre/home/harrell_lab/src/STAR-Fusion-WDL-hg38v22/` on the cluster (see Section 5), so
> the scripts are already available there — you do not need to clone them again. **If you
> are running this pipeline for a different project or lab, consider cloning the repository
> into a more neutral, shared location** and updating the repo paths in your commands
> accordingly.

### 8.4 Update `run_command.sh` and run

Point the run command at the **batch** WDL and batch inputs:

```bash
java -Dconfig.file=/lustre/home/harrell_lab/bulkRNASeq/config_wdl/COAD_PAAD_FusionAnalysis4Yang_260504/cromwell.conf \
  -jar ~/.cromwell/lib/cromwell-91.jar \
  run /lustre/home/harrell_lab/src/STAR-Fusion-WDL-hg38v22/star_fusion_hg38_batch_hpc_wf.wdl \
  --inputs batch_inputs.json 2>&1 | tee cromwell_run.log
```

Then launch it inside `tmux` exactly as in [Section 7.2](#72-run).

---

## 9. Monitoring

Multiple SLURM jobs will appear for batch runs (one per sample).

```bash
# SLURM job status
squeue -u $USER

# Job history and exit codes
sacct --user $USER --format=JobID,JobName,State,ExitCode,Elapsed,MaxRSS

# Watch the Cromwell log
tail -f cromwell_run.log
```

---

## 10. Verifying Outputs

On success, Cromwell prints `workflow finished with status 'Succeeded'`.

For batch runs, check the gathered `output_dir` (see [Section 8.2](#82-edit-batch_configjson)).
For single-sample runs, locate outputs within the execution tree:

```bash
find cromwell-executions/ -name "*.fusion_predictions.tsv*" 2>/dev/null
find cromwell-executions/ -name "*.coding_effect.tsv" 2>/dev/null
find cromwell-executions/ -name "*FusionInspector*abridged.tsv" 2>/dev/null
```

Key output files per sample:

| File | Description |
|------|-------------|
| `<sample>.star-fusion.fusion_predictions.tsv.gz` | Full fusion predictions |
| `<sample>.star-fusion.fusion_predictions.abridged.tsv.gz` | Abridged predictions |
| `star-fusion.fusion_predictions.abridged.coding_effect.tsv` | Coding effect (if enabled) |
| `finspector.FusionInspector.fusions.abridged.tsv` | FusionInspector results (if `validate`) |
| `<sample>.Chimeric.out.junction.gz` | STAR junction file (reusable for re-runs) |
| `<sample>.Log.final.out` | STAR alignment log |

---

## 11. Troubleshooting

When a run fails, inspect the Cromwell execution directory:

```bash
find cromwell-executions/ -name "docker_stderr" | xargs cat
find cromwell-executions/ -name "stdout"        | xargs cat
find cromwell-executions/ -name "rc"            | xargs cat   # return code
find cromwell-executions/ -name "script"        | xargs cat   # generated script
```

Common issues:

- **Exit code 127 (command not found):** Usually a bind-mount problem — the container
  can't see input files or the reference genome. Verify your `--bind` path covers every
  path in the inputs JSON; each location needs its own `--bind` entry.
- **Exit code 137 (OOM):** Increase `memory` in `cromwell.conf` or use a higher-memory
  partition. STAR-Fusion is memory-heavy during genome index loading.
- **"No such file or directory" for `cromwell.conf`:** Use the full absolute path with
  `-Dconfig.file=`.
- **`local_genome_dir` / genome not found inside container:** `local_genome_dir` is
  passed as a String (not a File), so Cromwell does not copy it — it must be reachable via
  the bind mount.
- **"File not found" for inputs:** Confirm the *real* path (`readlink -f`) is under a
  bound directory.
- **Import failure ("Could not find `star_fusion_workflow.wdl`"):** Cromwell resolves
  imports relative to the main WDL's location — point at the WDL inside its cloned repo
  directory, not a copy elsewhere.
- **Singularity image pull timeout:** The `trinityctat/starfusion:latest` image is several
  GB; the first pull is slow. If compute nodes lack internet, pre-pull on a node that has
  it:
  ```bash
  srun --partition=cpu-small --time=01:00:00 --pty /bin/bash
  module load singularity
  singularity pull docker://trinityctat/starfusion:latest
  ```
- **Disk space:** Check with `lfs quota -u $USER /lustre` and `df -h /lustre`.

---

## 12. Cleaning Up

Cromwell execution directories can be large. After verifying gathered results in your
output directory:

```bash
rm -rf cromwell-executions
```

---

## 13. Downstream Analysis and Filtering

The WDL pipeline produces **raw STAR-Fusion predictions only.** Turning those into
interpretable results requires downstream filtering and analysis — read-support and FFPM
thresholds, artifact and same-gene-family removal, priority annotation against known-fusion
databases, optional FusionInspector (Tier 2) validation, tumor fusion burden, cohort
enrichment, and expression analysis.

Scripts that perform this downstream work **for a specific COAD/CRC project** live in a
separate repository:

> **https://github.com/VCU-Bioinformatics-Core/fusion-gene-summary-report**

**Use these as a reference, not a turnkey step.** They document the process and the code is
reusable, but they were written for one dataset and its cohort definitions. Some steps
(specific filters, cohort configurations, hardcoded paths, priority tiers) may not apply to
your data. Review each step and adapt it to your project rather than running the scripts
blindly.

---

## 14. Running on Terra

The same pipeline imports into Terra using the **Terra/Cloud wrapper**
(`star_fusion_hg38_wf.wdl`). No cluster, Cromwell, or Singularity setup is required — Terra
manages execution.

1. Upload `star_fusion_hg38_wf.wdl` and `star_fusion_workflow.wdl` to your Terra workspace
   via the **Workflows** tab.
2. For a fresh run, populate `left_fq` / `right_fq` (or `fastq_pair_tar_gz`) in your data
   table.
3. For a junction re-run, populate `input_chimeric_junction` with the
   `Chimeric.out.junction.gz` path from a prior run and leave the FASTQ columns empty.

The Cloud wrapper takes `genome_plug_n_play_tar_gz` (default points to a GCS-hosted March
2021 build) and extracts the genome at runtime, so no pre-extraction is needed.

**Key differences from HPC:**

| Setting | Terra (Cloud) | Local HPC |
|---|---|---|
| WDL wrapper | `star_fusion_hg38_wf.wdl` | `star_fusion_hg38_hpc_wf.wdl` |
| Genome input | `genome_plug_n_play_tar_gz` (extracted at runtime) | `local_genome_dir` (pre-extracted directory) |
| Container runtime | Docker (managed) | Singularity via `module load` |
| Job scheduler | Google Cloud Pipelines API | SLURM |
| Reference / input paths | `gs://` or DRS URIs | Local (e.g., Lustre) paths |
| Preemptible / SSD / disk multipliers | Available | Not applicable |
| Submit method | API call | `sbatch` via submit script |

This guide does not provide step-by-step Terra submission instructions; the intent is to
document that the pipeline is Terra-compatible so it can be imported into future Terra
workspaces.

---

## 15. Credits and Citation

**Authors**

- **Brian Haas** (bhaas@broadinstitute.org) — Original STAR-Fusion pipeline
- **Amy Olex** (alolex@vcu.edu) — Junction-file input modifications, HPC execution support

This pipeline modifies the original STAR-Fusion Terra WDL from the
[Trinity CTAT project](https://github.com/STAR-Fusion/STAR-Fusion). If you use it, please
cite the original STAR-Fusion publication:

> Haas BJ, et al. *Accuracy assessment of fusion transcript detection via read-mapping and
> de novo fusion transcript assembly-based methods.* Genome Biology, 2019.

**AI disclosure:** Portions of the pipeline — including modified WDL code and documentation
— were developed with AI assistance. All AI-generated output was reviewed and validated
prior to use. Users are encouraged to validate pipeline behavior in their own environment
before production use.

---

## Cluster-Specific Notes (Apollo)

This workflow was developed and validated on the VCU Apollo HPC (SLURM, Lustre,
Singularity-CE 4.0.1). Adapting to another cluster typically involves: different partition
names (`sinfo -s`), different Singularity/Apptainer module names (`module avail`),
different base data paths for bind mounts, and localization method
(Apollo requires soft-link; some clusters support hard-link). On Apollo, `sbatch --wrap` is
not used — submit scripts are required, and this is already configured in the template.
