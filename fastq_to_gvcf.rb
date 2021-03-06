#!/usr/bin/env ruby
#
# analyze_sequence_to_snps.rb
# Created by Stuart Glenn on 2010-12-22
#
# == Synopsis
# Quick script to build the "analyze.sh" script used to run through the
# pipeline of steps for analysis going from fastq to vcf
#
# == Inputs
# The actual input to this script is just an job config file, which is used
# to easily specify all the parameters for many samples in on place. Future
# released might allow for single sample runs without the config file but for
# right now it is not a priority as we do analysis in bulk
#
# === Config File
# At this point the config file is a YAML hash of hashes. The hash is as follows
#
# {sample_id:  {is_paired: true/false, run: name, lane: number, :inputs: [fastq1,fastq2]}}
#
# It can also have a DEFAULT key with {bwa_ref: path_to_ref, gatk_ref: path_to_ref, snp_rod: path_to_rod}
# Those reference keys can be in each sample to override, if :snp_rod is missing we'll skip the covariate recalibration
#
#
# == Usage
#  analyze_sequence_to_snps.rb PATH_TO_CONFIG.YML PATH_TO_BASE_OUTPUT_DIR SAMPLE_ID [SAMPLE_ID]
#
# == Options
#  -h, --help             Display this help message
#  -v, --version          Display the version information
#  -V, --verbose          Increased verbosity of output
#  -c, --config FILE      Specify the configuration yaml file of options for analysis
#  -o, --output DIR       Specify the output directory prefix, all results will be saved under this directory
#
# ==Author
# Stuart Glenn <Stuart-Glenn@omrf.org>
#
# ==Copyright
#  Copyright (c) 2011-2017 Stuart Glenn, Oklahoma Medical Research Foundation. (OMRF)
#  All rights reserved.
#  Redistribution and use in source and binary forms, with or without
#  modification, are permitted provided that the following conditions are met:
#  1. Redistributions of source code must retain the above copyright
#     notice, this list of conditions and the following disclaimer.
#  2. Redistributions in binary form must reproduce the above copyright
#     notice, this list of conditions and the following disclaimer in the
#     documentation and/or other materials provided with the distribution.
#  3. All advertising materials mentioning features or use of this software
#     must display the following acknowledgement:
#     This product includes software developed by the OMRF
#  4. Neither the name of the Oklahoma Medical Research Foundation nor the
#     names of its contributors may be used to endorse or promote products
#     derived from this software without specific prior written permission.
#
#  THIS SOFTWARE IS PROVIDED BY COPYRIGHT HOLDERS AND CONTRIBUTORS ''AS IS'' AND ANY
#  EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
#  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
#  DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
#  DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
#  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
#  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
#  ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
#  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
#  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#

require 'yaml'
require 'erb'
require 'optparse'
require 'ostruct'
require 'json'


class Template

# Create instance
# *+default_config+ - Default options
# *+sample_name+ - Name of this sample
# *+data+ - Array of hashes of options for this sample, overrides anything in default_config, each is a library run
def initialize(default_config,sample_name,data)
  @default_config = default_config
  @sample_name = sample_name
  @data = data
end #initialize(default_config,sample_name,data)

def to_s
  ERB.new(script_template()).result(binding)
end

def scheduler_opts()
  @default_config[:opts][:scheduler_opts]
end

def haplocaller_opts(data)
  [data.first[:opts], @default_config[:opts]].each do |d|
    if d && d.has_key?(:haplotype_caller_opts)
      return d[:haplotype_caller_opts]
    end
  end
  return []
end

def tmp_dir_base_opt()
  base = @data.first[:opts][:tmp_dir_base] || @default_config[:opts][:tmp_dir_base] || nil
  if base
    "$(mktemp -d --suffix=.${SAMPLE}.$$ --tmpdir=\"#{base}\")"
  else
    "/tmp"
  end
end

def mode()
  self.class.which_mode(@data.first[:mode])
end

def self.which_mode(mode)
  case mode
  when /\Adna\z/i
    :dna
  when /\Arna\z/i
    :rna
  else
    raise "Unknown mode '#{data.first[:mode]}' for #{@sample_name}"
  end
end

def self.analysis_template(default_config,sample_name,data)
  case which_mode(data.first[:mode])
  when :dna
    DNAAnalysisTemplate.new(default_config,sample_name,data)
  when :rna
    RNAAnalysisTemplate.new(default_config,sample_name,data)
  end
end

def aligner_unload_load()
  case mode()
  when :dna
    "module unload bwa\nmodule load bwa/0.7.15\nmodule unload trimmomatic\nmodule load trimmomatic/0.35"
  when :rna
    "module unload star\nmodule load star/2.4.0h"
  end
end


def bash_header()
  <<-EOS
#!/bin/bash

cd "$( dirname "${BASH_SOURCE[0]}" )"

# be sure to start with a fresh & known enviroment (hopefully)
. /usr/local/analysis/lmod/lmod/init/bash
module unload samtools
module unload picard
module unload gatk
module unload fastqc
#{aligner_unload_load()}
module load samtools/1.3
module load picard/2.1.0
module load gatk/3.8-0-ge9d806836
module load fastqc/0.11.5

set -o pipefail
EOS
end
end

class UnalignedExtractTemplate < Template

def bam_to_fastq_for_run(run,bam_index)
  cmd = "bam2fastq --no-aligned --unaligned -o unaligned_fastq/#{@sample_name}_unaligned_#{run[:run]}_#{run[:lane]}\#_sequence.txt 03_first_bam/#{bam_index}.bam"
  cmd = "qsub #{scheduler_opts} -o logs -b y -j y -cwd -V -l h_vmem=8G -sync y -N a_#{@sample_name}_unaligned_#{bam_index} #{cmd}"
  return cmd
end

def script_template()
  <<-EOS
<%= bash_header() %>
module unload bam2fastq
module load bam2fastq/1.1.0

mkdir unaligned_fastq
<% @data.each_with_index do |run,i| %>
<%= bam_to_fastq_for_run(run,i) %>
if [ "$?" -ne "0" ]; then
echo "Failure extracting unalinged from bam <%= i %>"
exit 1
fi
<% end %>
qsub <%= scheduler_opts() %> -o logs -b y -V -j y -cwd -sync y -N a_<%=@sample_name%>_unaligned_gzip gzip -7 unaligned_fastq/*_sequence.txt

rm -rf 03_first_bam
  EOS
end
end

class AnalysisTemplate < Template

def fastq_file_list(sample_name,data)
  fastqs = []
  @fastq_shell_vars = {}
  @fastq_shell_vars_by_lane = []
  @input_sam_files = []
  letters = %w/A B C D E F G H I J K L M N O P Q R S T U V W X Y Z/
  data.each_with_index do |sequence,s_i|
    letter = letters[s_i]
    @fastq_shell_vars_by_lane << []
    sequence[:inputs].each_with_index do |input,i_i|
      prefix = "#{sample_name}_#{sequence[:run]}_#{sequence[:lane]}".downcase
      cleaned_prefix = "#{sample_name}_cleaned_#{sequence[:run]}_#{sequence[:lane]}".downcase
      pair_part = sequence[:is_paired] ? i_i+1 : 0
      shell_var = "FASTQ#{fastqs.size+1}"
      base_file = "#{cleaned_prefix}_#{pair_part}".downcase
      fastqs << "#{shell_var}=\"#{input}\""
      @fastq_shell_vars[shell_var] = {:path  => input, :paired => pair_part, :letter => letter, :base_file => base_file, :prefix => prefix}
      @fastq_shell_vars_by_lane[-1] << shell_var
      @input_sam_files << {:index => s_i, :b_index => pair_part}
    end
  end
  fastqs.join("\n")
end

def fastq_shell_vars()
  @fastq_shell_vars.keys.map{|v| "${#{v}}"}.join(" ")
end


def ordered_bam_inputs()
  @input_sam_files.map {|s| "#{s[:b_index]} ./00_inputs/#{s[:index]}.bam"}.join(" ")
end

def total_number_input_sequence_files
  @fastq_shell_vars.keys.size
end

def total_number_input_sequenced_lanes
  @fastq_shell_vars_by_lane.size
end

def link_fastq_inputs()
  #ln -s ${FASTQ1} 00_inputs/A_1.fastq
  lines = []
  @fastq_shell_vars.each do |var,data|
    lines << "ln -s ${#{var}} 00_inputs/#{data[:letter]}_#{data[:paired]}.fastq"
  end
  lines.join("\n")
end

def reference_for_data(data)
  data.first[:bwa_ref] || @default_config[:bwa_ref]
end

def alignment_command(sample_name,data)
end

def extract_unaligned(sample_name,data)

end

def input_sam_bam_files(prefix,suffix)
  cmd = ""
  # 02_bwa_alignment/A.sam 02_bwa_alignment/B.sam
  (0...@fastq_shell_vars_by_lane.size()).to_a.each do |i|
    cmd += "#{prefix}/#{i}.#{suffix} "
  end
  cmd
end

# output a -D SNP rod file, if we have a snp_rod
def opt_d_rod_path(data)
  if data.first[:snp_rod] || @default_config[:snp_rod]
    "-D ${GATK_DBSNP}"
  else
    ""
  end
end

# output a -L interval file, if we have an interval
def opt_l_interval(data)
  if data.first.has_key?(:interval_file)
    i = data.first[:interval_file]
    if (nil != i && "" != i && " " != i)
      return "-L #{i} -ip 50"
    end
  end
  return ""
end

def hs_metrics(sample_name,data)
  bam_dir = "."
  picard_mode,metric = if data.first.has_key?(:interval_file) then
                         [:hs,"CollectHsMetrics"]
                       else
                         [:wgs,"CollectWgsMetrics"]
                       end
  cmd="JAVA_MEM_OPTS=\"-Xmx24G\" sbatch #{scheduler_opts()} -N 1 -n 1 -c 1 --mem 36G \\
  -o logs/slurm-%x.%A.log -J a_#{sample_name}_alignment_summary \\\n"
  cmd+=<<EOS.chomp
  --wrap "picard #{metric} \\
  INPUT=#{bam_dir}/#{sample_name}.bam \\
  OUTPUT=#{bam_dir}/#{sample_name}_#{picard_mode}_metrics.txt \\
  VALIDATION_STRINGENCY=LENIENT \\
  REFERENCE_SEQUENCE=${GATK_REF}
EOS
  if :hs == picard_mode then
    cmd += " \\\n  BAIT_INTERVALS=#{data.first[:interval_file]} \\\n  TARGET_INTERVALS=#{data.first[:interval_file]}"
  end
  cmd+="\""
  return cmd
end

def variant_call(sample_name,data)
  if @default_config[:opts][:skip_gvcf]
    return ""
  end
  caller = ""
  bam_dir = "."
  if data.first.has_key?(:interval_file) then
    # if we have an interval file, just do it with that
    caller = ERB.new(<<-EOF
    # Finally Haplotypecaller in gVCF mode or is Gvcf mode

    export JAVA_MEM_OPTS="-Xmx24G"
    sbatch <%= scheduler_opts() %> -N 1 -n 1 -c 8 -o logs/slurm-%x.%A.log -W -J a_<%= sample_name %>_variants \\
    --mem 28 --wrap "gatk -T HaplotypeCaller \\
    -ERC GVCF -nct 8 -R ${GATK_REF} \\
    -I ./<%= bam_dir %>/<%= sample_name %>.bam -G StandardAnnotation -G AS_StandardAnnotation -G StandardHCAnnotation \\
    -o <%= sample_name %>.g.vcf.gz <%= opt_d_rod_path(data) %> <%= opt_l_interval(data) %> <%= haplocaller_opts(data).join(" ") %>"

    if [ "$?" -ne "0" ]; then
     echo "Failure GVCF"
     exit 1
    fi
EOF
    ).result(binding)
  else
    caller = gvcf_by_chr(sample_name,data)
  end

  return caller
end

def gvcf_by_chr(sample_name,data)
  chr_gvcfs = %w/1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y MT/.map do |chr|
    "-V 15_gvcf/by_chr/#{sample_name}-#{chr}.g.vcf.gz"
  end
  snprod = if data.first[:snp_rod] || @default_config[:snp_rod]
    "-d ${GATK_DBSNP}"
  else
    ""
  end
  extra_args = haplocaller_opts(data).map{|a| " -E \\\"#{a}\\\""}
  ERB.new(<<-EOF
  # GVCF by chr for better/faster/stronger
  mkdir 15_gvcf
  export JAVA_MEM_OPTS="-Xmx16G"
  qsub <%= scheduler_opts() %> -pe threaded 6 -t 1-25 -o logs -sync y -b y -V -j y -cwd -N a_<%= sample_name %>_gvcf_by_chr \\
  -l virtual_free=3G,mem_free=3G,h_vmem=7G haplocaller_qsub_tasked.rb -m 32 -r ${GATK_REF} <%= snprod %> \\
  -b 15_gvcf -p <%= sample_name %> -i ./<%= sample_name %>.bam <%= extra_args.join(" ") %>

  if [ "$?" -ne "0" ]; then
   echo "Failure GVCFing"
   exit 1
  fi

  export JAVA_MEM_OPTS="-Xmx24G"
  qsub <%= scheduler_opts() %> -o logs -sync y -b y -V -j y -cwd -N a_<%= sample_name %>_join_gvcf \\
  -l virtual_free=20G,mem_free=20G,h_vmem=30G gatk -T CombineGVCFs -R ${GATK_REF} \\
  <%= chr_gvcfs.join(" ") %> -o <%= sample_name %>.g.vcf.gz

  if [ "$?" -ne "0" ]; then
   echo "Failure joining gvcfs"
   exit 1
  fi
  rm -rf 15_gvcf

  EOF
  ).result(binding)
end

def known_indels_opts()
  if @default_config[:known_indels]
    return @default_config[:known_indels].map {|i| "-known #{i}"}.join(" ")
  else
    return ""
  end
end

def gatk_steps(sample_name,data)
end

def fastqc_bam(sample_name,data)
  if @default_config[:opts][:fastqc_bam]
    return <<-EOF
# fastqc info
mkdir qc
qsub #{scheduler_opts()} -p -1000 -l virtual_free=2G,h_vmem=4G -o logs -b y -V -j y -cwd -N a_#{sample_name}_qc fastqc -o qc ./#{sample_name}.bam
    EOF
  else
    return ""
  end
end

def fastqc_fastq(sample_name)
  if @default_config[:opts][:fastqc_fastq]
    return <<-EOF
# fastqc info
mkdir qc
qsub #{scheduler_opts()} -p -1000 -l virtual_free=2G,h_vmem=4G -o logs -b y -V -j y -cwd -N a_#{sample_name}_qc fastqc -o qc #{fastq_shell_vars()}
    EOF
  else
    return ""
  end
end

def indel_realign(sample_name,data)
  if @default_config[:opts][:skip_indel_realign]
    d = File.dirname(@current_input)
    @current_input = "07_realigned_bam/cleaned.bam"
    return <<-EOF
    ln -s #{d} 07_realigned_bam
    EOF
  else
    indel_realignment(sample_name,data)
  end
end

def indel_realgment_additional_opts
  ""
end

def bqsr_additional_opts
  indel_realgment_additional_opts() + " -SQQ 10 -SQQ 20 -SQQ 30 -SQQ 40"
end

def recalibrate_additional_opts()
  indel_realgment_additional_opts() + " --disable_indel_quals"
end

def indel_realignment(sample_name,data)
  compression = if data.first[:recalibration_known_sites] || @default_config[:recalibration_known_sites]
                  ""
                else
                  "--bam_compression 8"
                end
  f = @current_input
  @current_input = "07_realigned_bam/cleaned.bam"
  ERB.new(<<-EOF
    # Calculate intervals for realignment
    mkdir 06_intervals
    qsub <%= scheduler_opts() %> -pe threaded 6 -R y -o logs -sync y -b y -V -j y -cwd -N a_<%= sample_name %>_intervals \\
     -l virtual_free=1G,mem_free=1G,h_vmem=20G \\
     gatk -T RealignerTargetCreator -R ${GATK_REF} -I ./#{f} -o ./06_intervals/cleaned.intervals -nt 10 #{opt_l_interval(data)}

    if [ "$?" -ne "0" ]; then
     echo "Failure with target realigment creation"
     exit 1
    fi

    # Now realign & fix any mate info
    mkdir 07_realigned_bam
    unset JAVA_MEM_OPTS
    qsub <%= scheduler_opts() %> -o logs -sync y -b y -V -j y -cwd -N a_<%= sample_name %>_realign \\
     -l virtual_free=5G,mem_free=4G,h_vmem=8G \\
     gatk -T IndelRealigner <%= known_indels_opts() %> -R ${GATK_REF} -I ./04_dup_marked/cleaned.bam \\
     --targetIntervals ./06_intervals/cleaned.intervals -o ./07_realigned_bam/cleaned.bam --maxReadsInMemory 1000000 #{indel_realgment_additional_opts()} #{compression}

    if [ "$?" -ne "0" ]; then
     echo "Failure with indel realigmnent"
     exit 1
    fi
  EOF
  ).result(binding)
end

def mark_dupes_or_skip(sample_name,data)
  if @default_config[:opts][:skip_dupes]
    <<-EOF
sbatch #{scheduler_opts()} -t 0-12 -N 1 -n 1 -c 2 --mem 24 -o logs/slurm-%x.%A.log -W -J a_#{sample_name}_merge \\
 --wrap "picard MergeSamFiles #{input_sam_bam_files("INPUT=03_sorted_bams","bam")} OUTPUT=./04_dup_marked/cleaned.bam \\
 VALIDATION_STRINGENCY=LENIENT MAX_RECORDS_IN_RAM=6000000 CREATE_INDEX=True USE_THREADING=True COMPRESSION_LEVEL=8"

if [ "$?" -ne "0" ]; then
  echo "Failure with merging the sams"
  exit 1
fi
    EOF
  else
    mark_dupes(sample_name,data)
  end
end

def mark_dupes(sample_name,data)
  <<-EOF
sbatch #{scheduler_opts()} -t 0-12 -N 1 -n 1 -c 2 --mem 24 -o logs/slurm-%x.%A.log -W -J a_#{sample_name}_merge_mark_dups \\
  --wrap "picard MarkDuplicates #{input_sam_bam_files("INPUT=03_sorted_bams","bam")} \\
  OUTPUT=./04_dup_marked/cleaned.bam METRICS_FILE=./04_dup_marked/mark_dups_metrics.txt \\
  VALIDATION_STRINGENCY=LENIENT MAX_RECORDS_IN_RAM=6000000 CREATE_INDEX=True COMPRESSION_LEVEL=8"

if [ "$?" -ne "0" ]; then
  echo "Failure with marking the duplicates"
  exit 1
fi
  EOF
end

def covariate_or_final(sample_name,data)
  if data.first[:recalibration_known_sites] || @default_config[:recalibration_known_sites]
    covariate_recalibration(sample_name,data)
  else
    skip_covariate_recalibration(sample_name,data)
  end
end

def skip_covariate_recalibration(sample_name,data)
  <<-EOF
  mv ./07_realigned_bam/cleaned.bam ./#{sample_name}.bam
  mv ./07_realigned_bam/cleaned.bai ./#{sample_name}.bai
  EOF
end

def recalibration_known_sites()
  @default_config[:recalibration_known_sites].map {|s| "--knownSites #{s}"}.join(" ")
end

def covariate_recalibration(sample_name,data)
  ERB.new(<<-EOF
  # BaseRecalibrator
  mkdir 08_uncalibated_covariates
  unset JAVA_MEM_OPTS
  sbatch <%= scheduler_opts() %> -N 1 -n 1 -c 8 -o logs/slurm-%x.%A.log -W \\
    -J a_<%= sample_name %>_bqsr --mem 10 \\
    --wrap \"gatk -T BaseRecalibrator -R ${GATK_REF} <%= recalibration_known_sites() %> -I ./07_realigned_bam/cleaned.bam \\
    -o ./08_uncalibated_covariates/recal_data.grp -nct 10 <%= bqsr_additional_opts() %> <%= opt_l_interval(data) %>\"

  if [ "$?" -ne "0" ]; then
   echo "Failure counting covariates"
   exit 1
  fi

  mkdir 10_recalibrated_bam
  unset JAVA_MEM_OPTS
  sbatch <%= scheduler_opts() %> -N 1 -n 1 -c 8 -o logs/slurm-%x.%A.log -W \\
   -J a_<%= sample_name %>_recalibrate --mem 10 \\
   --wrap \"gatk -T PrintReads -R ${GATK_REF} -I ./07_realigned_bam/cleaned.bam -BQSR ./08_uncalibated_covariates/recal_data.grp \\
   -o ./10_recalibrated_bam/recalibrated.bam --bam_compression 8 -nct 10 <%= recalibrate_additional_opts %>\"

  if [ "$?" -ne "0" ]; then
   echo "Failure reclibrating bam"
   exit 1
  fi

  mv ./10_recalibrated_bam/recalibrated.bam ./<%= sample_name %>.bam
  mv ./10_recalibrated_bam/recalibrated.bai ./<%= sample_name %>.bam.bai
EOF
  ).result(binding)
end

# What strand call confidence for UnifiedGenotyper
# Will depend on if we are doing vqsr
def unified_genotyper_strand_call_conf(data)
  "30.0"
end

def script_template()
  return DATA.read
end

def path_variables()
  ERB.new(<<-EOF
GATK_BIN=$(which gatk)
GATK_BASE=$(dirname "${GATK_BIN}")"/.."
GATK_REF=<%= @data.first[:gatk_ref] || @default_config[:gatk_ref] %>
GATK_DBSNP=<%= @data.first[:snp_rod] || @default_config[:snp_rod] || 'nil' %>
EOF
).result(binding)
end

end

class DNAAnalysisTemplate < AnalysisTemplate

def gatk_steps(sample_name,data)
  @current_input = "04_dup_marked/cleaned.bam"
  <<-EOF
#{indel_realign(sample_name,data)}

#{covariate_or_final(sample_name,data)}
EOF
end

def reference_for_data(data)
  data.first[:bwa_ref] || @default_config[:bwa_ref]
end

def alignment_command(sample_name,data)
  cmd = "sbatch #{scheduler_opts()} -N 1 -c 12 --mem 48 -o logs/slurm-%x.%A.log -W -t 0-23 \\\n"
  cmd += " --array 1-#{total_number_input_sequenced_lanes()} -J a_#{sample_name}_bwa_alignment \\\n"
  cmd += " --wrap \"bwa-mem-slurm-array.rb -o 03_sorted_bams -r #{reference_for_data(data)}"
  cmd += " --source-env --download-timeout 1800 -t \"${TMP_DIR}\""
  if data.first[:opts].has_key?(:trimmomatic)
    data.first[:opts][:trimmomatic].each do |t|
      cmd += " --trim #{t}"
    end
  end
  @fastq_shell_vars_by_lane.each_with_index do |lane_shell_vars,index|
    if data[index][:is_paired]
      cmd += " paired"
    else
      cmd += " single"
    end

    # rg
    cmd += " '\"@RG\\\\tID:#{sample_name}_#{data[index][:run]}_s_#{data[index][:lane]}\\\\tSM:#{sample_name}\\\\tPL:Illumina\\\\tPU:#{data[index][:lane]}\"'"

    lane_shell_vars.each do |v|
      cmd += " \\\"${#{v}}\\\""
    end
  end
  cmd += "\""
  return cmd
end

end

class RNAAnalysisTemplate < AnalysisTemplate

def gatk_steps(sample_name,data)
  @current_input = "04_dup_marked/cleaned.bam"
  <<-EOF
#{split_trim_reassign_qualities(sample_name,data)}

#{indel_realign(sample_name,data)}

#{covariate_or_final(sample_name,data)}
EOF
end

def split_trim_reassign_qualities(sample_name,data)
  @current_input = "05_split-n-trim/cleaned.bam"
  <<-EOF
#Split'n'Trim for RNASeq
mkdir 05_split-n-trim

qsub #{scheduler_opts()} -R y -o logs -sync y -b y -V -j y -cwd -N a_#{sample_name}_splitntrim \\
 -l virtual_free=16G,mem_free=16G,h_vmem=20G \\
 gatk -T SplitNCigarReads -R ${GATK_REF} -I 04_dup_marked/cleaned.bam -o ./#{@current_input} \\
 -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS

 if [ "$?" -ne "0" ]; then
   echo "Failure with target realigment creation"
   exit 1
 fi
EOF
end

def path_variables()
  super() +
  ERB.new(<<-EOF
STAR_REF=<%= @data.first[:star_ref] || @default_config[:star_ref] %>
STAR_INDEX=<%= star_index() %>
EOF
).result(binding)
end

def star_index()
  @data.first[:star_index] || @default_config[:star_index]
end

def star_gtf()
  @data.first[:star_gtf] || @default_config[:star_gtf]
end

def alignment_command(sample_name,data)
  cmd = "qsub #{scheduler_opts()} -pe threaded 6 -l virtual_free=1G,mem_free=1G,h_vmem=48G -o logs -sync y \\\n"
  cmd += " -t 1-#{total_number_input_sequenced_lanes()} -b y -V -j y -cwd -N a_#{sample_name}_star_alignment \\\n"
  cmd += " star_qsub_tasked.rb -V -t \"${TMP_DIR}\" -o 03_sorted_bams -i ${STAR_INDEX} -r ${STAR_REF}"
  cmd += " -g #{star_gtf()}" if star_gtf()
  libs = []
  @fastq_shell_vars_by_lane.each_with_index do |lane_shell_vars,index|
    libs <<  {
      :paths => lane_shell_vars.map{|v| @fastq_shell_vars[v][:path]},
      :rg => {
      :id => "#{sample_name}_#{data[index][:run]}_s_#{data[index][:lane]}",
      :lb => "#{sample_name}_#{data[index][:run]}",
      :sm => sample_name,
      :pl => "Illumina",
      :pu => data[index][:lane],
    }}
  end
  j = libs.to_json.gsub(/"/,'\\"')
  cmd += " \\\n " + '\"' + "'" + j + "'" + '\"'
  return cmd
end

def indel_realgment_additional_opts
  "-U ALLOW_N_CIGAR_READS"
end

def bqsr_additional_opts
  indel_realgment_additional_opts()
end

def recalibrate_additional_opts
  indel_realgment_additional_opts()
end

end


class AnalysisTemplaterApp
VERSION       = "4.6.0"
REVISION_DATE = "20180111"
AUTHOR        = "Stuart Glenn <Stuart-Glenn@omrf.org>"
COPYRIGHT     = "Copyright (c) 2012-2017 Oklahoma Medical Research Foundation"

# Set up the app to roun
# *+args+ - command line ARGS formatted type array
# *+ios+ - optional hash to set stdin, stdout, & stderr
def initialize(args,ios = {})
  @stdin = ios[:stdin] || STDIN
  @stdout = ios[:stdout] || STDOUT
  @stderr = ios[:stderr] || STDERR

  @args = args
  set_default_options()
end #initialize(args,ios)

# Do the work
def run
  if !options_parsed?() || !options_valid?()
    @stderr.puts("")
    output_usage(@stderr)
    return 1
  end
  return run_real()
end #run

private

# Do the work of running this app
def run_real
  @default_config = (@config["DEFAULT"] || [
                     {:run => nil,
                      :lane => nil,
                      :bwa_ref => nil,
                      :star_ref => nil,
                      :star_gtf => nil,
                      :star_index => nil,
                      :gatk_ref => nil,
                      :snp_rod => nil,
                      :mode => nil,
                      :opts=>{
                        :skip_gvcf => false,
                        :skip_indel_realign => true,
                        :fastqc_bam => false,
                        :fastqc_fastq => false
                      }
                     }]).first

  @default_config[:opts].merge!(:scheduler_opts => @options.scheduler_opts) if '' != @options.scheduler_opts
  @default_config[:opts].merge!(:tmp_dir_base => @options.tmp_dir_base) if @options.tmp_dir_base
  statii = @options.samples.map do |sample_name|
    process_sample(sample_name)
  end
  return statii.max
end #run

# Make the dir, the analyze script and launch said script for sample_name
def process_sample(sample_name)
  data = @config[sample_name]
  data.each_with_index do |d,i|
    data[i] = @default_config.merge(d)
    case data[i][:mode]
    when /\ADNA\z/i
      if data[i][:bwa_ref] == nil
        @stderr.puts "Missing bwa reference for #{sample_name}"
        return 1
      end
    when /\Arna\z/i
      if data[i][:star_ref] == nil || data[i][:star_index] == nil
        @stderr.puts "Missing star reference for #{sample_name}"
        return 1
      end
    else
      @stderr.puts "Missing analysis mode for #{sample_name}"
      return 1
    end
  end
  if @options.debug
    puts data.inspect
    puts ""
    #exit 0
    puts Template.analysis_template(@default_config,sample_name,data)
    if (data.first.has_key?(:keep_unaligned) && data.first[:keep_unaligned]) then
      puts UnalignedExtractTemplate.new(@default_config,sample_name,data)
    end
    return 0
  end
  output_dir = File.join(@options.output_base,sample_name)

  unless Dir.mkdir(output_dir)
    @stderr.puts "Failed to make dir: #{output_dir} for #{sample_name}"
    return 1
  end

  script_file = File.join(output_dir,"analyze.sh")
  File.open(script_file,"w") do |f|
    f.puts Template.analysis_template(@default_config,sample_name,data)
  end

  if (data.first.has_key?(:keep_unaligned) && data.first[:keep_unaligned]) then
    extract_script_file = File.join(output_dir,"extract_unaligned.sh")
    File.open(extract_script_file,"w") do |f|
      f.puts UnalignedExtractTemplate.new(@default_config,sample_name,data)
    end
  end

  return_dir = Dir.pwd
  unless Dir.chdir(output_dir)
    @stderr.puts "Failed to change to dir: #{output_dir} for #{sample_name}"
    return 1
  end

  Dir.mkdir("logs")

  # We sleep a random amount to avoid overloading SGE with a billion jobs right away
  sleep(rand(@options.delay))
  cmd = %W(sbatch) + @options.scheduler_opts.split(/ /) + %W(-o logs/slurm-%x.%A.log -W - t 0-23 -N 1 -c 1 --mem 4 -J a_#{sample_name}_full ./analyze.sh)
  cmd = %w(./analyze.sh) if @options.run_local
  @stdout.puts(cmd.join(" "))
  system(*cmd)
  status = $?.exitstatus

  Dir.chdir(return_dir)
  return status
end #process_sample(sample_name)

# Parse the command line options, returning false on errors
def options_parsed?
  opts = OptionParser.new() do |o|
    o.on('-v','--version') { output_version($stdout); exit(0) }
    o.on('-h','--help') { output_help($stdout); exit(0) }
    o.on('-V', '--verbose')    { @options.verbose = true }
    o.on('-D', '--debug')    { @options.debug = true }
    o.on('-l', '--local')    { @options.run_local = true }

    o.on("-d","--delay", "=REQUIRED") do |amount|
      @options.delay = amount.to_i
    end

    o.on("-c","--config", "=REQUIRED") do |conf_file|
      @options.config_file = conf_file
    end

    o.on("-o","--output", "=REQUIRED") do |output_destination|
      @options.output_base = output_destination
    end

    o.on("-s","--scheduler", "=REQUIRED") do |qopts|
      @options.scheduler_opts = qopts
   end

    o.on("-t","--tmp", "=REQUIRED") do |topts|
      @options.tmp_dir_base = topts
    end
  end

  opts.parse!(@args) rescue return false
  @options.samples = @args
  return true
end #options_parsed?

# Open/parase/load the YAML config file
def loaded_config?
  begin
    @config = YAML::load(File.open(@options.config_file))
  rescue => err
    @stderr.puts "Error loading config file: #{err}"
    return false
  end
  return true
end #loaded_config?

# Test that the given output base dir exists & is writeabale
def output_base_valid?
  if File.exists?(@options.output_base)
    if !File.directory?(@options.output_base)
      @stderr.puts("Final output location, #{@options.output_base}, exists and is not a directory")
      return false
    elsif !File.writable?(@options.output_base)
      @stderr.puts("Final output location, #{@options.output_base}, is not writable")
      return false
    end
  else
    Dir.mkdir(@options.output_base)
  end
  return true
end #output_base_valid?

# Did the user provide us with additional args to use as sample ids
def have_sample_ids?
  if @options.samples && @options.samples.size > 0
    return true
  end
  @stderr.puts "Missing sample id(s) to processor"
  return false
end #have_sample_ids?

# Are all the ids they gave in the config file actually
def sample_ids_in_config?
  # check them all first
  @options.samples.each do |s|
    unless @config[s]
      @stderr.puts "Can't find sample '#{s}' in config"
      return false
    end
  end
  return true
end #sample_ids_in_config?

# Makes sure the options we have make sense
def options_valid?
  loaded_config? &&
  output_base_valid? &&
  have_sample_ids? &&
  sample_ids_in_config?
end #options_valid?

# Set our default options
def set_default_options()
  @options = OpenStruct.new(
    :output_base => nil,
    :verbose => false,
    :config_file  => nil,
    :samples => nil,
    :delay => 30,
    :scheduler_opts => ''
  )
end #set_default_options()

# Just print the usage message
def output_usage(out)
  out.puts <<-EOF
analyze_sequence_to_snps.rb -c PATH_TO_CONFIG.YML -o PATH_TO_BASE_OUTPUT_DIR SAMPLE_ID [SAMPLE_ID]

Options
-h, --help             Display this help message
-v, --version          Display the version information
-V, --verbose          Increased verbosity of output
-d, --delay INT        Delay a random between 0 & INT before submitting job. Default 30 seonds
-c, --config FILE      Specify the configuration yaml file of options for analysis
-o, --output DIR       Specify the output directory prefix, all results will be saved under this directory
-l, --local            Run the analyze script locally, not with initial SGE submit
-q, --qsub OPTS        Additional options given to each qsub call
  EOF
end

def output_version(out)
  out.puts "#{File.basename(__FILE__)} Version: #{VERSION} Released: #{REVISION_DATE}"
end

def output_help(out)
  output_version(out)
  out.puts ""
  output_usage(out)
end

end

if $0 == __FILE__
exit(AnalysisTemplaterApp.new(ARGV.clone).run())
end

__END__
<%=
bash_header()
%>

# It is easier to use variables for these paths
<%= path_variables() %>

SAMPLE="<%= @sample_name %>"

<%=
fastq_file_list(@sample_name,@data)
%>

function final_clean() {
rm -rf "$TMP_DIR"
}
TMP_DIR=<%= tmp_dir_base_opt() %>
if [ "$?" -ne "0" ]; then
  echo "Failed to make our temp work dir"
  exit 1
fi
trap final_clean EXIT
export GATK_JAVA_OPTS="-Djava.io.tmpdir=${TMP_DIR}"

if [ ! -e 04_dup_marked/cleaned.bam ]; then
<%= fastqc_fastq(@sample_name) %>

export JAVA_MEM_OPTS="-Xmx16G"

# Now take those & actually generate the alignment SAM output, paired or single
mkdir 03_sorted_bams
<%=
  alignment_command(@sample_name,@data)
%>

if [ "$?" -ne "0" ]; then
  echo "Failure with alignment"
  exit 1
fi

# Now we might have had many input SAMs, so let us merge those all into a single BAM using picard
# While we do that we shall also sort it, mark possible duplicates & make an index
mkdir 04_dup_marked
<%= mark_dupes_or_skip(@sample_name,@data) %>

<%= (@data.first.has_key?(:keep_unaligned) && @data.first[:keep_unaligned]) ? "": "rm -rf 03_sorted_bams" %>

if [ "$PRE_GATK_ONLY" == "Y" ]; then
  rm -rf 00_inputs \
  01_bwa_aln_sai \
  02_bwa_alignment <%= (@data.first.has_key?(:keep_unaligned) && @data.first[:keep_unaligned]) ? "": "03_sorted_bams" %> \
  "${TMP_DIR}"

  rm -f qc/*.zip

  touch finished.txt
  exit 0
fi
else
  rm -f finished.txt
  TMP_DIR=<%= tmp_dir_base_opt() %>
  if [ "$?" -ne "0" ]; then
    echo "Failed to make our temp work dir"
    exit 1
  fi
fi #if 04_dup_marked/cleaned.bam already existed
# start here if pre_gatk already done

<%= gatk_steps(@sample_name,@data) %>

<%= fastqc_bam(@sample_name,@data) %>

# Clean up after ourselves
rm -rf 00_inputs \
01_bwa_aln_sai \
02_bwa_alignment <%= (@data.first.has_key?(:keep_unaligned) && @data.first[:keep_unaligned]) ? "": "03_sorted_bams" %> \
04_dup_marked \
05_split-n-trim \
06_intervals \
07_realigned_bam \
08_uncalibated_covariates \
10_recalibrated_bam

# Get some summary stats on the alignment off in the background
<%= hs_metrics(@sample_name,@data) %>

<%= variant_call(@sample_name,@data) %>

<%=
  if (@data.first.has_key?(:keep_unaligned) && @data.first[:keep_unaligned]) then
    cmd = "qsub #{scheduler_opts()} -o logs -b y -V -j y -cwd -m e -N a_#{@sample_name}_unaligned_extract ./extract_unaligned.sh"
    cmd
  end
%>

rm -f qc/*.zip
rm -rf "${TMP_DIR}"

touch finished.txt
