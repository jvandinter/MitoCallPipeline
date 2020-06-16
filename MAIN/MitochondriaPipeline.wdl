## Version 1.2
## 1.0: Working order
## 1.1: Optimised runtime for faster running
##      All output to one folder
## 1.2: Added lineage-based blacklist

import "AlignAndCall.wdl" as AlignAndCall

  workflow MitochondriaPipeline {

  meta {
    description: "Takes in fully aligned hg38 bam and outputs VCF of SNP/Indel calls on the mitochondria."
  }


  ## Input names
  File wgs_aligned_input_bam_or_cram
  File wgs_aligned_input_bam_or_cram_index
  String sample_name
  String contig_name = "MT"
  Float autosomal_coverage

  # Read length used for optimization only. If this is too small CollectWgsMetrics might fail,
  # but the results are not affected by this number. Default = 151.
  Int? max_read_length

  # Full reference is only required if starting with a CRAM (BAM doesn't need these files)
  File? hg38_fasta
  File? hg38_index
  File? hg38_dict

  ## Reference genome files
  File mt_dict
  File mt_fasta
  File mt_fasta_index
  File mt_amb
  File mt_ann
  File mt_bwt
  File mt_pac
  File mt_sa
  File blacklisted_sites
  File blacklisted_sites_M 
  File blacklisted_sites_M_i
  File blacklisted_sites_L 
  File blacklisted_sites_L_i
  File blacklisted_sites_N 
  File blacklisted_sites_N_i

  #Shifted reference is used for calling the control region (edge of mitochondria reference).
  #This solves the problem that BWA doesn't support alignment to circular contigs.
  File mt_shifted_dict
  File mt_shifted_fasta
  File mt_shifted_fasta_index
  File mt_shifted_amb
  File mt_shifted_ann
  File mt_shifted_bwt
  File mt_shifted_pac
  File mt_shifted_sa

  File shift_back_chain

  File control_region_shifted_reference_interval_list
  File non_control_region_interval_list

  String? m2_extra_args
  String? m2_filter_extra_args
  Float? vaf_filter_threshold
  Float? f_score_beta

  ## Common resource claims

  # Versions of required tools
  String module_picard_version
  String module_gatk_version
  String module_samtools_version
  String module_r_version

  Boolean compress_output_vcf = false
  #String? tmp_mon_or_default = select_first([tmp_mon, "& pid=$!;trap \"kill $pid 2> /dev/null\" EXIT;while kill -0 $pid 2> /dev/null; do df -lh |grep $TMPDIR >>$HOME/tmpUse/$JOB_ID;sleep 60;done"])

  parameter_meta {
    wgs_aligned_input_bam_or_cram: "Full WGS hg38 bam or cram"
    autosomal_coverage: "Median coverage of full input bam"
    out_vcf: "Final VCF of mitochondrial SNPs and INDELs"
    vaf_filter_threshold: "Hard threshold for filtering low VAF sites"
    f_score_beta: "F-Score beta balances the filtering strategy between recall and precision. The relative weight of recall to precision."
    contig_name: "Name of mitochondria contig in reference that wgs_aligned_input_bam_or_cram is aligned to"
    sample_name: "Name of file in final output vcf"
  }

  call SubsetBamToChrM {
    input:
      input_bam = wgs_aligned_input_bam_or_cram,
      input_bai = wgs_aligned_input_bam_or_cram_index,
      contig_name = contig_name,
      hg38_fasta = hg38_fasta,
      hg38_index = hg38_index,
      hg38_dict = hg38_dict,
      module_gatk_version = module_gatk_version,
  }

  call RevertSam {
    input:
      input_bam = SubsetBamToChrM.output_bam,
      module_picard_version = module_picard_version,
      module_samtools_version = module_samtools_version,
  }

  call AlignAndCall.AlignAndCall as AlignAndCall {
    input:
      unmapped_bam = RevertSam.unmapped_bam,
      autosomal_coverage = autosomal_coverage,
      sample_name = sample_name,
      mt_dict = mt_dict,
      mt_fasta = mt_fasta,
      mt_fasta_index = mt_fasta_index,
      mt_amb = mt_amb,
      mt_ann = mt_ann,
      mt_bwt = mt_bwt,
      mt_pac = mt_pac,
      mt_sa = mt_sa,
      mt_shifted_dict = mt_shifted_dict,
      mt_shifted_fasta = mt_shifted_fasta,
      mt_shifted_fasta_index = mt_shifted_fasta_index,
      mt_shifted_amb = mt_shifted_amb,
      mt_shifted_ann = mt_shifted_ann,
      mt_shifted_bwt = mt_shifted_bwt,
      mt_shifted_pac = mt_shifted_pac,
      mt_shifted_sa = mt_shifted_sa,
      shift_back_chain = shift_back_chain,
      m2_extra_args = m2_extra_args,
      m2_filter_extra_args = m2_filter_extra_args,
      vaf_filter_threshold = vaf_filter_threshold,
      f_score_beta = f_score_beta,
      compress_output_vcf = compress_output_vcf,
      max_read_length = max_read_length,
      module_gatk_version = module_gatk_version,
      module_picard_version = module_picard_version,
      module_r_version = module_r_version,
  }

  # This is a temporary task to handle "joint calling" until Mutect2 can produce a GVCF.
  # This provides coverage at each base so low coverage sites can be considered ./. rather than 0/0.
  call CoverageAtEveryBase {
    input:
      input_bam_regular_ref = AlignAndCall.mt_aligned_bam,
      input_bam_regular_ref_index = AlignAndCall.mt_aligned_bai,
      input_bam_shifted_ref = AlignAndCall.mt_aligned_shifted_bam,
      input_bam_shifted_ref_index = AlignAndCall.mt_aligned_shifted_bai,
      shift_back_chain = shift_back_chain,
      control_region_shifted_reference_interval_list = control_region_shifted_reference_interval_list,
      non_control_region_interval_list = non_control_region_interval_list,
      ref_fasta = mt_fasta,
      ref_fasta_index = mt_fasta_index,
      ref_dict = mt_dict,
      shifted_ref_fasta = mt_shifted_fasta,
      shifted_ref_fasta_index = mt_shifted_fasta_index,
      shifted_ref_dict = mt_shifted_dict,
      module_picard_version = module_picard_version,
      module_r_version = module_r_version,
  }

  output {
    File subset_bam = SubsetBamToChrM.output_bam
    File subset_bai = SubsetBamToChrM.output_bai
    File mt_aligned_bam = AlignAndCall.mt_aligned_bam
    File mt_aligned_bai = AlignAndCall.mt_aligned_bai
    File out_vcf = AlignAndCall.out_vcf
    File out_vcf_index = AlignAndCall.out_vcf_index
    File duplicate_metrics = AlignAndCall.duplicate_metrics
    File coverage_metrics = AlignAndCall.coverage_metrics
    File theoretical_sensitivity_metrics = AlignAndCall.theoretical_sensitivity_metrics
    File contamination_metrics = AlignAndCall.contamination_metrics
    File base_level_coverage_metrics = CoverageAtEveryBase.table
    Int mean_coverage = AlignAndCall.mean_coverage
    String major_haplogroup = AlignAndCall.major_haplogroup
    Float contamination = AlignAndCall.contamination
  }
}

task SubsetBamToChrM {
  File input_bam
  File input_bai
  String contig_name
  String basename = basename(basename(input_bam, ".cram"), ".bam")
  String module_gatk_version
  File? hg38_fasta
  File? hg38_index
  File? hg38_dict

  meta {
    description: "Subsets a whole genome bam to just Mitochondria reads"
  }
  parameter_meta {
    ref_fasta: "Reference is only required for cram input. If it is provided ref_fasta_index and ref_dict are also required."
    }

  command <<<
    set -e
    module load ${module_gatk_version}

    java -Dpicard.useLegacyParser=false -Djava.io.tmpdir=$TMPDIR -jar $GATK4 \
    PrintReads \
      ${"-R " + hg38_fasta} \
      -L ${contig_name} \
      --read-filter MateOnSameContigOrNoMappedMateReadFilter \
      --read-filter MateUnmappedAndUnmappedReadFilter \
      -I ${input_bam} \
      -O ${basename}.bam
  >>>

  runtime {
    cpu : "1"
    memory : "4 GB"
    tmpspace_gb : "4"
    wallclock : "1:00:00"
  }
  output {
    File output_bam = "${basename}.bam"
    File output_bai = "${basename}.bai"
  }
}
task RevertSam {
  File input_bam
  String basename = basename(input_bam, ".bam")
  String module_picard_version
  String module_samtools_version

  meta {
    description: "Removes alignment information while retaining recalibrated base qualities and original alignment tags"
  }
  parameter_meta {
    input_bam: "aligned bam"
  }

  command <<<
    module load ${module_picard_version}
    module load ${module_samtools_version}

    java -Dpicard.useLegacyParser=false -Djava.io.tmpdir=$TMPDIR -Xms4G -Xmx12G -jar $PICARD \
    RevertSam \
    --INPUT ${input_bam} \
    --OUTPUT_BY_READGROUP false \
    --OUTPUT ${basename}.bam \
    --VALIDATION_STRINGENCY LENIENT \
    --ATTRIBUTE_TO_CLEAR FT \
    --ATTRIBUTE_TO_CLEAR CO \
    --SORT_ORDER queryname \
    --RESTORE_ORIGINAL_QUALITIES false
  >>>

  runtime {
    cpu: "1"
    memory: "16 GB"
    wallclock: "1:00:00"
    tmpspace_gb: "10"
  }
  output {
    File unmapped_bam = "${basename}.bam"
  }
}

task CoverageAtEveryBase {
  File input_bam_regular_ref
  File input_bam_regular_ref_index
  File input_bam_shifted_ref
  File input_bam_shifted_ref_index
  File shift_back_chain
  File control_region_shifted_reference_interval_list
  File non_control_region_interval_list
  File ref_fasta
  File ref_fasta_index
  File ref_dict
  File shifted_ref_fasta
  File shifted_ref_fasta_index
  File shifted_ref_dict
  String module_picard_version
  String module_r_version

  meta {
    description: "Remove this hack once there's a GVCF solution."
  }

  command <<<
    set -e
    module load ${module_picard_version}
    module load ${module_r_version}

    java -Dpicard.useLegacyParser=false -Djava.io.tmpdir=$TMPDIR -jar $PICARD \
    CollectHsMetrics \
      -I ${input_bam_regular_ref} \
      -R ${ref_fasta} \
      --PER_BASE_COVERAGE non_control_region.tsv \
      -O non_control_region.metrics \
      -TI ${non_control_region_interval_list} \
      -BI ${non_control_region_interval_list} \
      --COVERAGE_CAP 20000 \
      --SAMPLE_SIZE 1 

    java -Dpicard.useLegacyParser=false -Djava.io.tmpdir=$TMPDIR -jar $PICARD \
    CollectHsMetrics \
      -I ${input_bam_shifted_ref} \
      -R ${shifted_ref_fasta} \
      --PER_BASE_COVERAGE control_region_shifted.tsv \
      -O control_region_shifted.metrics \
      -TI ${control_region_shifted_reference_interval_list} \
      -BI ${control_region_shifted_reference_interval_list} \
      --COVERAGE_CAP 20000 \
      --SAMPLE_SIZE 1

    R --vanilla <<CODE
      shift_back = function(x) {
        if (x < 8570) {
          return(x + 8000)
        } else {
          return (x - 8569)
        }
      }

      control_region_shifted = read.table("control_region_shifted.tsv", header=T)
      shifted_back = sapply(control_region_shifted[,"pos"], shift_back)
      control_region_shifted[,"pos"] = shifted_back

      beginning = subset(control_region_shifted, control_region_shifted[,'pos']<8000)
      end = subset(control_region_shifted, control_region_shifted[,'pos']>8000)

      non_control_region = read.table("non_control_region.tsv", header=T)
      combined_table = rbind(beginning, non_control_region, end)
      write.table(combined_table, "per_base_coverage.tsv", row.names=F, col.names=T, quote=F, sep="\t")

    CODE
  >>>

  runtime {
    cpu: "1"
    memory: "4 GB"
    wallclock: "1:00:00"
    tmpspace_gb: "10"
  }
  output {
    File table = "per_base_coverage.tsv"
  }
}