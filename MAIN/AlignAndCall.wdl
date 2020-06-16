## Version 1.2
## 1.0: Working order
## 1.1: Optimised runtime for faster running
##      All output to one folder
## 1.2: Added lineage-based blacklist

import "AlignmentPipeline.wdl" as AlignAndMarkDuplicates

  workflow AlignAndCall {

  meta {
    description: "Takes in unmapped bam and outputs VCF of SNP/Indel calls on the mitochondria."
  }


  File unmapped_bam
  Float? autosomal_coverage
  String sample_name

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

  # Shifted reference is used for calling the control region (edge of mitochondria reference).
  # This solves the problem that BWA doesn't support alignment to circular contigs.
  File mt_shifted_dict
  File mt_shifted_fasta
  File mt_shifted_fasta_index
  File mt_shifted_amb
  File mt_shifted_ann
  File mt_shifted_bwt
  File mt_shifted_pac
  File mt_shifted_sa

  # Module versions
  String module_gatk_version
  String module_python_version
  String module_picard_version
  String module_bwa_version
  String module_samtools_version
  String module_haplochecker_location
  String module_r_version

  File shift_back_chain

  String? m2_extra_args
  String? m2_filter_extra_args
  Float? vaf_filter_threshold
  Float? f_score_beta
  Boolean compress_output_vcf

  # Read length used for optimization only. If this is too small CollectWgsMetrics might fail, but the results are not
  # affected by this number. Default is 151.
  Int? max_read_length

  parameter_meta {
    unmapped_bam: "Unmapped and subset bam, optionally with original alignment (OA) tag"
  }

  call AlignAndMarkDuplicates.AlignmentPipeline as AlignToMt {
    input:
      input_bam = unmapped_bam,
      mt_dict = mt_dict,
      mt_fasta = mt_fasta,
      mt_fasta_index = mt_fasta_index,
      mt_amb = mt_amb,
      mt_ann = mt_ann,
      mt_bwt = mt_bwt,
      mt_pac = mt_pac,
      mt_sa = mt_sa,
      module_bwa_version = module_bwa_version,
      module_picard_version = module_picard_version,
  }

  call AlignAndMarkDuplicates.AlignmentPipeline as AlignToShiftedMt {
    input:
      input_bam = unmapped_bam,
      mt_dict = mt_shifted_dict,
      mt_fasta = mt_shifted_fasta,
      mt_fasta_index = mt_shifted_fasta_index,
      mt_amb = mt_shifted_amb,
      mt_ann = mt_shifted_ann,
      mt_bwt = mt_shifted_bwt,
      mt_pac = mt_shifted_pac,
      mt_sa = mt_shifted_sa,
      module_bwa_version = module_bwa_version,
      module_picard_version = module_picard_version,
  }

  call CollectWgsMetrics {
    input:
      input_bam = AlignToMt.mt_aligned_bam,
      input_bam_index = AlignToMt.mt_aligned_bai,
      ref_fasta = mt_fasta,
      ref_fasta_index = mt_fasta_index,
      read_length = max_read_length,
      coverage_cap = 100000,
      module_picard_version = module_picard_version,
      module_r_version = module_r_version
  }

  call GetContamination {
    input:
      input_bam = AlignToMt.mt_aligned_bam,
      input_bam_index = AlignToMt.mt_aligned_bai,
      ref_fasta = mt_fasta,
      ref_fasta_index = mt_fasta_index,
      module_python_version = module_python_version,
      module_haplochecker_location = module_haplochecker_location,
  }


  call M2 as CallMt {
    input:
      input_bam = AlignToMt.mt_aligned_bam,
      input_bai = AlignToMt.mt_aligned_bai,
      ref_fasta = mt_fasta,
      ref_fai = mt_fasta_index,
      ref_dict = mt_dict,
      compress = compress_output_vcf,
      module_gatk_version = module_gatk_version,
      # Everything is called except the control region.
      m2_extra_args = select_first([m2_extra_args, ""]) + " -L chrM:576-16024 ",
  }

  call M2 as CallShiftedMt {
    input:
      input_bam = AlignToShiftedMt.mt_aligned_bam,
      input_bai = AlignToShiftedMt.mt_aligned_bai,
      ref_fasta = mt_shifted_fasta,
      ref_fai = mt_shifted_fasta_index,
      ref_dict = mt_shifted_dict,
      compress = compress_output_vcf,
      module_gatk_version = module_gatk_version,
      # Everything is called except the control region.
      m2_extra_args = select_first([m2_extra_args, ""]) + " -L chrM:8025-9144 ",
  }

  call LiftoverAndCombineVcfs {
    input:
      shifted_vcf = CallShiftedMt.raw_vcf,
      vcf = CallMt.raw_vcf,
      ref_fasta = mt_fasta,
      ref_fasta_index = mt_fasta_index,
      ref_dict = mt_dict,
      shift_back_chain = shift_back_chain,
      module_picard_version = module_picard_version,
  }

  call MergeStats {
    input:
      shifted_stats = CallShiftedMt.stats,
      non_shifted_stats = CallMt.stats,
      module_gatk_version = module_gatk_version,
  }
  call Haplolist {
    input:
      haplogroup = GetContamination.major_hg,
      blacklisted_sites = blacklisted_sites,
      blacklisted_sites_L = blacklisted_sites_L,
      blacklisted_sites_L_i = blacklisted_sites_L_i,
      blacklisted_sites_M = blacklisted_sites_M,
      blacklisted_sites_M_i = blacklisted_sites_M_i,
      blacklisted_sites_N = blacklisted_sites_N,
      blacklisted_sites_N_i = blacklisted_sites_N_i,
      module_python_version = module_python_version
  }

  call Filter {
    input:
      raw_vcf = LiftoverAndCombineVcfs.final_vcf,
      raw_vcf_index = LiftoverAndCombineVcfs.final_vcf_index,
      raw_vcf_stats = MergeStats.stats,
      sample_name = sample_name,
      ref_fasta = mt_fasta,
      ref_fai = mt_fasta_index,
      ref_dict = mt_dict,
      compress = compress_output_vcf,
      m2_extra_filtering_args = m2_filter_extra_args,
      max_alt_allele_count = 4,
      contamination = GetContamination.minor_level,
      autosomal_coverage = autosomal_coverage,
      vaf_filter_threshold = vaf_filter_threshold,
      used_blacklist = Haplolist.used_blacklist,
      used_index = Haplolist.used_index,
      f_score_beta = f_score_beta,
      module_gatk_version = module_gatk_version,
  }

  output {
    File mt_aligned_bam = AlignToMt.mt_aligned_bam
    File mt_aligned_bai = AlignToMt.mt_aligned_bai
    File mt_aligned_shifted_bam = AlignToShiftedMt.mt_aligned_bam
    File mt_aligned_shifted_bai = AlignToShiftedMt.mt_aligned_bai
    File out_vcf = Filter.filtered_vcf
    File out_vcf_index = Filter.filtered_vcf_idx
    File duplicate_metrics = AlignToMt.duplicate_metrics
    File coverage_metrics = CollectWgsMetrics.metrics
    File theoretical_sensitivity_metrics = CollectWgsMetrics.theoretical_sensitivity
    File contamination_metrics = GetContamination.contamination_file
    Int mean_coverage = CollectWgsMetrics.mean_coverage
    String major_haplogroup = GetContamination.major_hg
    Float contamination = GetContamination.minor_level
  }
}

task GetContamination {
  File input_bam
  File input_bam_index
  File ref_fasta
  File ref_fasta_index
  String module_haplochecker_location
  String module_python_version
  Int qual = 20
  Int map_qual = 30
  Float vaf = 0.01

  String basename = basename(input_bam, ".bam")

  meta {
    description: "Uses Haplochecker to estimate levels of contamination in mitochondria"
  }
  parameter_meta {
    input_bam: "Bam aligned to chrM"
    ref_fasta: "chrM reference"
  }
  command <<<
    set -e

    java -jar ${module_haplochecker_location} \
      haplochecker \
      --in ${input_bam} \
      --ref ${ref_fasta} \
      --out haplochecker_out \
      --QUAL ${qual} \
      --MAPQ ${map_qual} \
      --VAF ${vaf} 

module load ${module_python_version}

python3 <<CODE

import csv

with open("haplochecker_out/${basename}.contamination.txt") as output:
    reader = csv.DictReader(output, delimiter='\t')
    for row in reader:
        print(row["MajorHG"], file=open("major_hg.txt", 'w'))
        print(row["MajorLevel"], file=open("major_level.txt", 'w'))
        print(row["MinorHG"], file=open("minor_hg.txt", 'w'))
        print(row["MinorLevel"], file=open("minor_level.txt", 'w'))
CODE
  >>>
  runtime {
    cpu: "1"
    memory: "4 GB"
    tmpspace_gb: "10"
    wallclock: "1:00:00"
  }
  output {
    File contamination_file = "haplochecker_out/${basename}.contamination.txt"
    String major_hg = read_string("major_hg.txt")
    Float major_level = read_float("major_level.txt")
    String minor_hg = read_string("minor_hg.txt")
    Float minor_level = read_float("minor_level.txt")
  }
}

task CollectWgsMetrics {
  String module_picard_version
  String module_r_version
  File input_bam
  File input_bam_index
  File ref_fasta
  File ref_fasta_index
  Int? read_length
  Int? coverage_cap

  Int read_length_for_optimization = select_first([read_length, 151])
  Float ref_size = size(ref_fasta, "GB") + size(ref_fasta_index, "GB")
  Int disk_size = ceil(size(input_bam, "GB") + ref_size) + 20

  meta {
    description: "Collect coverage metrics"
  }
  parameter_meta {
    read_length: "Read length used for optimization only. If this is too small CollectWgsMetrics might fail. Default is 151."
  }

  command <<<
    set -e
    module load ${module_picard_version}
    module load ${module_r_version}

    java -Dpicard.useLegacyParser=false -Djava.io.tmpdir=$TMPDIR -jar $PICARD \
      CollectWgsMetrics \
      --INPUT ${input_bam} \
      --VALIDATION_STRINGENCY SILENT \
      --REFERENCE_SEQUENCE ${ref_fasta} \
      --OUTPUT metrics.txt \
      --USE_FAST_ALGORITHM true \
      --READ_LENGTH ${read_length_for_optimization} \
      --${"COVERAGE_CAP=" + coverage_cap} \
      --INCLUDE_BQ_HISTOGRAM true \
      --THEORETICAL_SENSITIVITY_OUTPUT theoretical_sensitivity.txt 

    R --vanilla <<CODE
      df = read.table("metrics.txt",skip=6,header=TRUE,stringsAsFactors=FALSE,sep='\t',nrows=1)
      write.table(floor(df[,"MEAN_COVERAGE"]), "mean_coverage.txt", quote=F, col.names=F, row.names=F)
    CODE
  >>>
  runtime {
    cpu: "1"
    memory: "4 GB"
    tmpspace_gb : "10"
    wallclock : "1:00:00"
  }
  output {
    File metrics = "metrics.txt"
    File theoretical_sensitivity = "theoretical_sensitivity.txt"
    Int mean_coverage = read_int("mean_coverage.txt")
  }
}

task LiftoverAndCombineVcfs {
  File shifted_vcf
  File vcf
  String basename = basename(shifted_vcf, ".vcf")

  File ref_fasta
  File ref_fasta_index
  File ref_dict

  File shift_back_chain
  String module_picard_version

  # runtime

  Float ref_size = size(ref_fasta, "GB") + size(ref_fasta_index, "GB")
  Int disk_size = ceil(size(shifted_vcf, "GB") + ref_size) + 20

  meta {
    description: "Lifts over shifted vcf of control region and combines it with the rest of the chrM calls."
  }
  parameter_meta {
    shifted_vcf: "VCF of control region on shifted reference"
    vcf: "VCF of the rest of chrM on original reference"
    ref_fasta: "Original (not shifted) chrM reference"
    shift_back_chain: "Chain file to lift over from shifted reference to original chrM"
  }
  command<<<
    set -e

    module load ${module_picard_version}
    java -Dpicard.useLegacyParser=false -Djava.io.tmpdir=$TMPDIR -jar $PICARD \
      LiftoverVcf \
      -I ${shifted_vcf} \
      -O ${basename}.shifted_back.vcf \
      -R ${ref_fasta} \
      --CHAIN ${shift_back_chain} \
      --REJECT ${basename}.rejected.vcf 

    java -Dpicard.useLegacyParser=false -Djava.io.tmpdir=$TMPDIR -jar $PICARD \
      MergeVcfs \
      -I ${basename}.shifted_back.vcf \
      -I ${vcf} \
      -O ${basename}.final.vcf
    >>>
    runtime {
     cpu : "1"
     memory : "4 GB"
     tmpspace_gb : "10"
     wallclock : "1:00:00"
    }
    output{
        # rejected_vcf should always be empty
        File rejected_vcf = "${basename}.rejected.vcf"
        File final_vcf = "${basename}.final.vcf"
        File final_vcf_index = "${basename}.final.vcf.idx"
    }
}

task M2 {
  File ref_fasta
  File ref_fai
  File ref_dict
  File input_bam
  File input_bai
  String module_gatk_version
  String? m2_extra_args
  Boolean? make_bamout
  Boolean compress
  File? gga_vcf
  File? gga_vcf_idx
  # runtime

  String output_vcf = "raw" + if compress then ".vcf.gz" else ".vcf"
  String output_vcf_index = output_vcf + if compress then ".tbi" else ".idx"
  #Float ref_size = size(ref_fasta, "GB") + size(ref_fai, "GB")
  #Int disk_size = ceil(size(input_bam, "GB") + ref_size) + 20

  # Mem is in units of GB but our command and memory runtime values are in MB
  # Not required -JVD
  #Int machine_mem = if defined(mem) then mem * 1000 else 3500
  #Int command_mem = machine_mem - 500

  meta {
    description: "Mutect2 for calling Snps and Indels"
  }
  parameter_meta {
    input_bam: "Aligned Bam"
    gga_vcf: "VCF for genotype given alleles mode"
  }
  command <<<
      set -e

      module load ${module_gatk_version}
      # We need to create these files regardless, even if they stay empty
      touch bamout.bam

      java -Djava.io.tmpdir=$TMPDIR -Xmx8G -jar $GATK4 \
        Mutect2 \
        -R ${ref_fasta} \
        -I ${input_bam} \
        ${"--genotyping-mode GENOTYPE_GIVEN_ALLELES --alleles " + gga_vcf} \
        -O ${output_vcf} \
        ${true='--bam-output bamout.bam' false='' make_bamout} \
        ${m2_extra_args} \
        --annotation StrandBiasBySample \
        --mitochondria-mode \
        --max-reads-per-alignment-start 75 \
        --max-mnp-distance 0
  >>>
  runtime {
      cpu: "2"
      memory: "12 GB"
      tmpspace_gb : "10"
      wallclock : "1:00:00"
  }
  output {
      File raw_vcf = "${output_vcf}"
      File raw_vcf_idx = "${output_vcf_index}"
      File stats = "${output_vcf}.stats"
      File output_bamOut = "bamout.bam"
  }
}

task Haplolist {
  File blacklisted_sites
  File blacklisted_sites_M 
  File blacklisted_sites_M_i
  File blacklisted_sites_L 
  File blacklisted_sites_L_i
  File blacklisted_sites_N 
  File blacklisted_sites_N_i
  String module_python_version
  String haplogroup

  command <<<
    set -e
    module load ${module_python_version}


    python3 <<CODE

    import csv
    haplo_L = ['L']
    haplo_M = ['C','D','E','G','M','Q','Z']
    haplo_N = ['A','B','F','H','I','J','K','N','O','P','R','S','T','U','V','W','X','Y']
    haploletter = "${haplogroup}"[0]

    def search_haplo (hl, hg):
        for i in hg:
            if i == hl:
                return 0
        return 1
    with open('blacklist_used.txt', 'w') as output:
        with open('index_used.txt', 'w') as index:
            if search_haplo(haploletter, haplo_L) == 0:
                output.write('${blacklisted_sites_L}')
                index.write('${blacklisted_sites_L_i}') 
            elif search_haplo(haploletter, haplo_M) == 0:
                output.write('${blacklisted_sites_M}')
                index.write('${blacklisted_sites_M_i}')
            elif search_haplo(haploletter, haplo_N) == 0:
                output.write('${blacklisted_sites_N}')
                index.write('${blacklisted_sites_N_i}')
            else: 
                output.write('${blacklisted_sites}')
                index.write('${blacklisted_sites}' + '.tbi')
    CODE
  >>>
  
  runtime {
    cpu: "1"
    memory: "4 GB"
    tmpspace_gb: "1"
    wallclock: "1:00:00"
  }
  output {
    File used_blacklist = read_string("blacklist_used.txt")
	File used_index = read_string("index_used.txt")
  }
}

task Filter {
  File ref_fasta
  File ref_fai
  File ref_dict
  File raw_vcf
  File raw_vcf_index
  File raw_vcf_stats
  Boolean compress
  Float? vaf_cutoff
  String sample_name

  String? m2_extra_filtering_args
  Int max_alt_allele_count
  Float contamination
  Float? autosomal_coverage
  Float? vaf_filter_threshold
  Float? f_score_beta

  File used_blacklist
  File used_index

  String module_gatk_version

  # runtime

  String output_vcf = sample_name + if compress then ".vcf.gz" else ".vcf"
  String output_vcf_index = output_vcf + if compress then ".tbi" else ".idx"

  meta {
    description: "Mutect2 Filtering for calling Snps and Indels"
  }
  parameter_meta {
      autosomal_coverage: "Median coverage of the autosomes for filtering potential polymorphic NuMT variants"
      vaf_filter_threshold: "Hard cutoff for minimum allele fraction. All sites with VAF less than this cutoff will be filtered."
      f_score_beta: "F-Score beta balances the filtering strategy between recall and precision. The relative weight of recall to precision."
      haplogroup: "Determines the blacklist with common variants per mitochondrial haplogroup to use."
  }
  command <<<
      set -e

       module load ${module_gatk_version}

      # We need to create these files regardless, even if they stay empty
      touch bamout.bam

        java -Xmx10G -Djava.io.tmpdir=$TMPDIR -Xms6g -Xss4096k -jar $GATK4 \
          FilterMutectCalls \
          -V ${raw_vcf} \
          -R ${ref_fasta} \
          -O filtered.vcf \
          --stats ${raw_vcf_stats} \
          ${m2_extra_filtering_args} \
          --max-alt-allele-count ${max_alt_allele_count} \
          --mitochondria-mode \
          ${"--autosomal-coverage " + autosomal_coverage} \
          ${"--min-allele-fraction " + vaf_filter_threshold} \
          ${"--f-score-beta " + f_score_beta} \
          --contamination-estimate ${contamination}

        java -Djava.io.tmpdir=$TMPDIR -jar $GATK4 \
          VariantFiltration \
          -V filtered.vcf \
          -O ${output_vcf} \
          --mask ${used_blacklist} \
          --read-index ${used_index} \
          --mask-name "blacklisted_site"

  >>>
  runtime {
      cpu: "1"
      memory: "14 GB"
      tmpspace_gb : "1"
      wallclock : "1:00:00"
  }
  output {
      File filtered_vcf = "${output_vcf}"
      File filtered_vcf_idx = "${output_vcf_index}"
  }
}

task MergeStats {
  File shifted_stats
  File non_shifted_stats
  String module_gatk_version

  command{
    set -e

    module load ${module_gatk_version}

    java -Djava.io.tmpdir=$TMPDIR -jar $GATK4 \
    MergeMutectStats \
    --stats ${shifted_stats} \
    --stats ${non_shifted_stats} \
    -O raw.combined.stats
  }
  output {
    File stats = "raw.combined.stats"
  }
  runtime {
      memory: "4 GB"
      tmpspace_gb: "1"
      wallclock: "1:00:00"
  }
}
