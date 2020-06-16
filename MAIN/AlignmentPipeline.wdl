## Version 1.2
## 1.0: Working order
## 1.1: Optimised runtime for faster running
##      All output to one folder
## 1.2: Added lineage-based blacklist

workflow AlignmentPipeline {

  meta {
    description: "Uses BWA to align unmapped bam and marks duplicates."
  }

  File input_bam
  String basename = basename(input_bam, ".bam")

  # FASTA files
  File mt_dict
  File mt_fasta
  File mt_fasta_index
  File mt_amb
  File mt_ann
  File mt_bwt
  File mt_pac
  File mt_sa

  #Module versions
  String module_bwa_version
  String module_picard_version

  ## Common resource claims

  parameter_meta {
    input_bam: "Input is an unaligned subset bam of NuMT and chrM reads and their mates. All reads must be paired."
    mt_aligned_bam: "Output is aligned duplicate marked coordinate sorted bam."
  }

  call AlignAndMarkDuplicates {
    input:
      input_bam = input_bam,
      module_bwa_version = module_bwa_version,
      module_picard_version = module_picard_version,
      output_bam_basename = basename + ".realigned",
      ref_fasta = mt_fasta,
      ref_fasta_index = mt_fasta_index,
      ref_dict = mt_dict,
      ref_amb = mt_amb,
      ref_ann = mt_ann,
      ref_bwt = mt_bwt,
      ref_pac = mt_pac,
      ref_sa = mt_sa,
  }
  output {
    File mt_aligned_bam = AlignAndMarkDuplicates.output_bam
    File mt_aligned_bai = AlignAndMarkDuplicates.output_bam_index
    File duplicate_metrics = AlignAndMarkDuplicates.duplicate_metrics
  }
}
task AlignAndMarkDuplicates {
  File input_bam
  String bwa_commandline = "bwa mem -K 100000000 -p -v 3 -t 2 -Y $bash_ref_fasta"
  String module_bwa_version
  String module_picard_version
  String output_bam_basename
  File ref_fasta
  File ref_fasta_index
  File ref_dict
  File ref_amb
  File ref_ann
  File ref_bwt
  File ref_pac
  File ref_sa
  String basename = basename(input_bam, ".bam")
  String metrics_filename = basename + ".metrics"

  meta {
    description: "Aligns with BWA and MergeBamAlignment, then Marks Duplicates. Outputs a coordinate sorted bam."
  }
  parameter_meta {
    input_bam: "Unmapped bam"
    module_bwa_version: "BWA version to be added to header of aligned bam"
  }
  command <<<
    module load ${module_picard_version}
    module load ${module_bwa_version}

    set -o pipefail
    set -e

    bash_ref_fasta=${ref_fasta}
    java -Dpicard.useLegacyParser=false -jar $PICARD \
      SamToFastq \
      --INPUT ${input_bam} \
      --FASTQ /dev/stdout \
      --INTERLEAVE true \
      --NON_PF true | \
    ${bwa_commandline} /dev/stdin -  2> >(tee ${output_bam_basename}.bwa.stderr.log >&2) | \
    java -Dpicard.useLegacyParser=false -Djava.io.tmpdir=$TMPDIR -jar $PICARD \
      MergeBamAlignment \
      -VALIDATION_STRINGENCY SILENT \
      -EXPECTED_ORIENTATIONS FR \
      -ATTRIBUTES_TO_RETAIN X0 \
      -ATTRIBUTES_TO_REMOVE NM \
      -ATTRIBUTES_TO_REMOVE MD \
      --ALIGNED_BAM /dev/stdin \
      --UNMAPPED_BAM ${input_bam} \
      --OUTPUT mba.bam \
      --REFERENCE_SEQUENCE ${ref_fasta} \
      --PAIRED_RUN true \
      --SORT_ORDER "unsorted" \
      --IS_BISULFITE_SEQUENCE false \
      --ALIGNED_READS_ONLY false \
      --CLIP_ADAPTERS false \
      --MAX_RECORDS_IN_RAM 2000000 \
      --ADD_MATE_CIGAR true \
      --MAX_INSERTIONS_OR_DELETIONS -1 \
      --PRIMARY_ALIGNMENT_STRATEGY MostDistant \
      --PROGRAM_RECORD_ID "bwamem" \
      --PROGRAM_GROUP_VERSION "${module_bwa_version}" \
      --PROGRAM_GROUP_COMMAND_LINE "${bwa_commandline}" \
      --PROGRAM_GROUP_NAME "bwamem" \
      --UNMAPPED_READ_STRATEGY COPY_TO_TAG \
      --ALIGNER_PROPER_PAIR_FLAGS true \
      --UNMAP_CONTAMINANT_READS true \
      --ADD_PG_TAG_TO_READS false 

    java -Dpicard.useLegacyParser=false -Djava.io.tmpdir=$TMPDIR -jar $PICARD \
      MarkDuplicates \
      --INPUT mba.bam \
      --OUTPUT md.bam \
      --METRICS_FILE ${metrics_filename} \
      --VALIDATION_STRINGENCY SILENT \
      --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \
      --ASSUME_SORT_ORDER "queryname" \
      --CLEAR_DT "false" \
      --ADD_PG_TAG_TO_READS false 

      java -Dpicard.useLegacyParser=false -Djava.io.tmpdir=$TMPDIR -jar $PICARD \
        SortSam \
        --INPUT md.bam \
        --OUTPUT ${output_bam_basename}.bam \
        --SORT_ORDER "coordinate" \
        --CREATE_INDEX true \
        --MAX_RECORDS_IN_RAM 300000
  >>>
  runtime {
    cpu: "1"
    memory: "16 GB"
    tmpspace_gb: "10"
    wallclock: "1:00:00"
  }
  output {
    File output_bam = "${output_bam_basename}.bam"
    File output_bam_index = "${output_bam_basename}.bai"
    File bwa_stderr_log = "${output_bam_basename}.bwa.stderr.log"
    File duplicate_metrics = "${metrics_filename}"
  }
}

