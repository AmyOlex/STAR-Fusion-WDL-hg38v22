version 1.0

workflow star_fusion_workflow {
  input {

    String sample_id

    File genome_plug_n_play_tar_gz

    # input data options - FASTQs
    File? left_fq
    File? right_fq
    File? fastq_pair_tar_gz

    # input data options - Chimeric junction file (alternative to FASTQs)
    File? input_chimeric_junction  # Chimeric.out.junction or Chimeric.out.junction.gz

    # STAR-Fusion parameters
    String? fusion_inspector  # inspect or validate
    Boolean examine_coding_effect = false
    Boolean coord_sort_bam = false
    Float min_FFPM = 0.1

    # runtime params
    String docker = "trinityctat/starfusion:latest"
    Int num_cpu = 12
    Float fastq_disk_space_multiplier = 3.25
    String memory = "50G"
    Float genome_disk_space_multiplier = 2.5
    Int preemptible = 2
    Float extra_disk_space = 10
    Boolean use_ssd = true

  }

  call star_fusion {
    input:
      left_fq                      = left_fq,
      right_fq                     = right_fq,
      fastq_pair_tar_gz            = fastq_pair_tar_gz,
      input_chimeric_junction      = input_chimeric_junction,
      genome                       = genome_plug_n_play_tar_gz,
      sample_id                    = sample_id,
      examine_coding_effect        = examine_coding_effect,
      coord_sort_bam               = coord_sort_bam,
      min_FFPM                     = min_FFPM,

      preemptible                  = preemptible,
      docker                       = docker,
      cpu                          = num_cpu,
      memory                       = memory,
      extra_disk_space             = extra_disk_space,
      fastq_disk_space_multiplier  = fastq_disk_space_multiplier,
      genome_disk_space_multiplier = genome_disk_space_multiplier,
      fusion_inspector             = fusion_inspector,
      use_ssd                      = use_ssd
  }

  output {
    File fusion_predictions           = star_fusion.fusion_predictions
    File fusion_predictions_abridged  = star_fusion.fusion_predictions_abridged
    File? junction                    = star_fusion.junction
    File? bam                         = star_fusion.bam
    File? sj                          = star_fusion.sj

    File? coding_effect               = star_fusion.coding_effect
    Array[File]? extract_fusion_reads = star_fusion.extract_fusion_reads

    File? star_log_final              = star_fusion.star_log_final

    File? fusion_inspector_validate_fusions_abridged = star_fusion.fusion_inspector_validate_fusions_abridged
    File? fusion_inspector_validate_web              = star_fusion.fusion_inspector_validate_web

    File? fusion_inspector_inspect_fusions_abridged  = star_fusion.fusion_inspector_inspect_fusions_abridged
    File? fusion_inspector_inspect_web               = star_fusion.fusion_inspector_inspect_web
  }
}


task star_fusion {
  input {
    String sample_id

    File? left_fq
    File? right_fq
    File? fastq_pair_tar_gz

    File? input_chimeric_junction

    File genome

    String? fusion_inspector
    Boolean examine_coding_effect
    Boolean coord_sort_bam
    Float min_FFPM

    Int preemptible
    String docker
    Int cpu
    String memory
    Float extra_disk_space
    Float fastq_disk_space_multiplier
    Float genome_disk_space_multiplier
    Boolean use_ssd
  }

  Float junction_size = if defined(input_chimeric_junction) then size(input_chimeric_junction, "GB") else 0
  Float fastq_size    = (if defined(left_fq) then size(left_fq, "GB") else 0) +
                        (if defined(right_fq) then size(right_fq, "GB") else 0) +
                        (if defined(fastq_pair_tar_gz) then size(fastq_pair_tar_gz, "GB") else 0)
  Float input_disk    = if defined(input_chimeric_junction)
                        then junction_size * 2
                        else fastq_disk_space_multiplier * fastq_size

  command <<<
    set -ex
    shopt -s nullglob

    mkdir -p ~{sample_id}

    # -----------------------------------------------------------------------
    # Build input params depending on mode
    # -----------------------------------------------------------------------
    if [[ ! -z "~{input_chimeric_junction}" ]]; then

      # --- Junction file input mode ---
      chimeric_junc="~{input_chimeric_junction}"
      if [[ "${chimeric_junc}" == *.gz ]]; then
        gunzip -c "${chimeric_junc}" > chimeric.out.junction
        chimeric_junc="chimeric.out.junction"
      fi
      read_params="-J ${chimeric_junc}"

    else

      # --- FASTQ input mode (original behavior) ---
      if [[ ! -z "~{fastq_pair_tar_gz}" ]]; then
        # untar the fq pair
        mv ~{fastq_pair_tar_gz} reads.tar.gz
        tar xvf reads.tar.gz

        left_fq=(*_1.fastq* *_1.fq*)

        if [[ ! -z "*_2.fastq*" ]] || [[ ! -z "*_2.fq*" ]]; then
          right_fq=(*_2.fastq* *_2.fq*)
        fi
      else
        left_fq="~{left_fq}"
        right_fq="~{right_fq}"
      fi

      if [[ -z "${left_fq[0]}" && -z "${right_fq[0]}" ]]; then
        echo "Error, not finding fastq files here"
        ls -ltr
        exit 1
      fi

      left_fqs=$(IFS=, ; echo "${left_fq[*]}")

      read_params="--left_fq ${left_fqs}"
      if [[ "${right_fq[0]}" != "" ]]; then
        right_fqs=$(IFS=, ; echo "${right_fq[*]}")
        read_params="${read_params} --right_fq ${right_fqs}"
      fi

    fi

    # -----------------------------------------------------------------------
    # Unpack genome
    # -----------------------------------------------------------------------
    mkdir -p genome_dir
    tar xf ~{genome} -C genome_dir --strip-components 1

    # -----------------------------------------------------------------------
    # Run STAR-Fusion
    # -----------------------------------------------------------------------
    /usr/local/src/STAR-Fusion/STAR-Fusion \
      --genome_lib_dir `pwd`/genome_dir/ctat_genome_lib_build_dir \
      ${read_params} \
      --output_dir ~{sample_id} \
      --CPU ~{cpu} \
      ~{"--FusionInspector " + fusion_inspector} \
      ~{true='--examine_coding_effect' false='' examine_coding_effect} \
      ~{"--min_FFPM " + min_FFPM}

    # -----------------------------------------------------------------------
    # Sort / rename output BAM (only produced in FASTQ mode)
    # -----------------------------------------------------------------------
    if [ -f ~{sample_id}/Aligned.out.bam ]; then
      if [ "~{coord_sort_bam}" == "true" ]; then
        samtools sort -@~{cpu} -o ~{sample_id}.STAR.aligned.coordsorted.bam ~{sample_id}/Aligned.out.bam
        samtools index ~{sample_id}.STAR.aligned.coordsorted.bam
      else
        mv ~{sample_id}/Aligned.out.bam ~{sample_id}.STAR.aligned.UNsorted.bam
      fi
    fi

    # -----------------------------------------------------------------------
    # Rename / compress remaining outputs
    # -----------------------------------------------------------------------
    mv ~{sample_id}/star-fusion.fusion_predictions.tsv ~{sample_id}.star-fusion.fusion_predictions.tsv \
      && gzip ~{sample_id}.star-fusion.fusion_predictions.tsv
    mv ~{sample_id}/star-fusion.fusion_predictions.abridged.tsv ~{sample_id}.star-fusion.fusion_predictions.abridged.tsv \
      && gzip ~{sample_id}.star-fusion.fusion_predictions.abridged.tsv

    if [ -f ~{sample_id}/Chimeric.out.junction ]; then
      mv ~{sample_id}/Chimeric.out.junction ~{sample_id}.Chimeric.out.junction && gzip ~{sample_id}.Chimeric.out.junction
    fi
    if [ -f ~{sample_id}/SJ.out.tab ]; then
      mv ~{sample_id}/SJ.out.tab ~{sample_id}.SJ.out.tab && gzip ~{sample_id}.SJ.out.tab
    fi
    if [ -f ~{sample_id}/Log.final.out ]; then
      mv ~{sample_id}/Log.final.out ~{sample_id}.Log.final.out
    fi

    gzip -c ~{sample_id}/star-fusion.preliminary/star-fusion.fusion_candidates.preliminary \
      > ~{sample_id}.star-fusion.fusion_candidates.preliminary.tsv.gz

  >>>

  output {
    File  fusion_predictions             = "~{sample_id}.star-fusion.fusion_predictions.tsv.gz"
    File  fusion_predictions_abridged    = "~{sample_id}.star-fusion.fusion_predictions.abridged.tsv.gz"
    File  preliminary_fusion_predictions = "~{sample_id}.star-fusion.fusion_candidates.preliminary.tsv.gz"

    File? junction       = "~{sample_id}.Chimeric.out.junction.gz"
    File? bam            = if (coord_sort_bam) then "~{sample_id}.STAR.aligned.coordsorted.bam" else "~{sample_id}.STAR.aligned.UNsorted.bam"
    File? sj             = "~{sample_id}.SJ.out.tab.gz"
    File? star_log_final = "~{sample_id}.Log.final.out"

    File? coding_effect              = "~{sample_id}/star-fusion.fusion_predictions.abridged.coding_effect.tsv"
    Array[File] extract_fusion_reads = glob("~{sample_id}/star-fusion.fusion_evidence_*.fq")

    File? fusion_inspector_validate_fusions_abridged = "~{sample_id}/FusionInspector-validate/finspector.FusionInspector.fusions.abridged.tsv"
    File? fusion_inspector_validate_web              = "~{sample_id}/FusionInspector-validate/finspector.fusion_inspector_web.html"
    File? fusion_inspector_inspect_fusions_abridged  = "~{sample_id}/FusionInspector-inspect/finspector.FusionInspector.fusions.abridged.tsv"
    File? fusion_inspector_inspect_web               = "~{sample_id}/FusionInspector-inspect/finspector.fusion_inspector_web.html"
  }

  runtime {
    preemptible: preemptible
    disks: "local-disk " + ceil(input_disk + size(genome, "GB") * genome_disk_space_multiplier + extra_disk_space) + " " + (if use_ssd then "SSD" else "HDD")
    docker: docker
    cpu: cpu
    memory: memory
  }
}
