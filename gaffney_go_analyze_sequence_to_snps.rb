#!/usr/bin/env ruby1.9

require 'yaml'
require 'erb'

def fastq_file_list(sample_name,data)
  fastqs = []
  @fastq_shell_vars = {}
  @fastq_shell_vars_by_lane = []
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
      path = "`pwd`\"/#{prefix}/#{base_file}.fastq\""
      fastqs << "#{shell_var}=#{path}"
      @fastq_shell_vars[shell_var] = {:path  => path, :paired => pair_part, :letter => letter, :base_file => base_file, :prefix => prefix}
      @fastq_shell_vars_by_lane[-1] << shell_var
    end
  end
  fastqs.join("\n")
end

def fastq_shell_vars()
  @fastq_shell_vars.keys.map{|v| "${#{v}}"}.join(" ")
end

def ordered_fastq_inputs()
  line = ""
  @fastq_shell_vars_by_lane.flatten.each do |input|
    line += "${#{input}} "
  end
  line
end

def gzip_original_fastq(sample_name)
  #qsub -o logs -b y -V -j y -cwd -q all.q -N <%= sample_name %>_gzip_1 gzip --fast ${FASTQ1}
  cmds = []
  @fastq_shell_vars_by_lane.flatten.each_with_index do |input,i|
    cmds << "qsub -o logs -b y -V -j y -cwd -q all.q -N #{sample_name}_gzip_#{i} gzip --fast ${#{input}}"
    cmds << "qsub -o logs -b y -V -j y -cwd -q all.q -N #{sample_name}_gzip_#{i}_rejects gzip --fast #{@fastq_shell_vars[input][:prefix]}/rejects.txt" if 0==i%2
  end
  cmds.join("\n")
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

def bwa_aligment_command(sample_name,data)
  cmd = "qsub -o logs -sync y -t 1-#{total_number_input_sequenced_lanes()} -b y -V -j y -cwd -q all.q -N #{sample_name}_bwa_alignment bwa_sampese_qsub_tasked.rb 02_bwa_alignment /Volumes/hts_core/Shared/homo_sapiens_36.1/chr_fixed/bwa_indexed/hg18.fa"
  @fastq_shell_vars_by_lane.each_with_index do |lane_shell_vars,index|
    if data[index][:is_paired]
      cmd += " paired"
    else
      cmd += " single"
    end

    # rg
    cmd += " '\"@RG\\\\tID:#{sample_name}_#{data[index][:run]}_s_#{data[index][:lane]}\\\\tSM:#{sample_name}\\\\tPL:Illumina\\\\tPU:#{data[index][:lane]}\"'"

    # 01_bwa_aln_sai file(s)
    lane_shell_vars.each do |v|
      cmd += " 01_bwa_aln_sai/#{@fastq_shell_vars[v][:base_file]}.sai"
    end
    # fastq file(s)
    lane_shell_vars.each do |v|
      cmd += " ${#{v}}"
    end
  end
  return cmd  
end

def default_rg(sample_name,data)
  return "" if data.first[:is_paired]
  data = data.first
  "--default_read_group #{sample_name}_#{data[:run]}_s_#{data[:lane]} --default_platform Illumina"
end

def input_sam_bam_files(prefix,suffix)
  cmd = ""
  # 02_bwa_alignment/A.sam 02_bwa_alignment/B.sam
  (0...@fastq_shell_vars_by_lane.size()).to_a.each do |i|
    cmd += "#{prefix}/#{i}.#{suffix} "
  end
  cmd
end

def clean_commands(sample_name,data)
  #clean_sample.rb -r ${RUN} -l ${LANE} -s <%= sample_name %> -b . [--single] INPUT_SEQUENCE
  cleans = []
  data.each_with_index do |sequence,s_i|
    cmd = "clean_sample.rb -s #{sample_name} -r #{sequence[:run]} -l #{sequence[:lane]} -b ."
    unless sequence[:is_paired]
      cmd += " --single-end"
    end
    cmd += " #{sequence[:inputs].join(" ").gsub(/\\/,"\\\\\\")}"
    cleans << "qsub -o logs -sync y -b y -V -j y -cwd -q all.q -N #{sample_name}_clean_#{s_i+1} #{cmd}"
  end
  cleans.join("\n")
end

output_base = ARGV.shift
config = YAML::load(File.open(ARGV.shift))

samples = ARGV.clone

script_template = <<-EOF
#!/bin/bash

module load bwa/0.5.9rc1
module load samtools/0.1.12
module load picard/1.36
module load gatk/1.0.4705
module load fastqc/0.7.2
module load btangs/1.2.0

GATK_REF=/Volumes/hts_core/Shared/homo_sapiens_36.1/chr_fixed/hg18.fasta
GATK_DBSNP=/Volumes/hts_core/Shared/dbsnp/dbsnp_129_hg18.rod

GATK_BIN=`which gatk`
GATK_BASE=`dirname ${GATK_BIN}`"/.."

SAMPLE="<%= sample_name %>"

<%=
  fastq_file_list(sample_name,data)
%>

# clean
<%=
  clean_commands(sample_name,data)
%>

if [ "$?" -ne "0" ]; then
  echo -e "Failure with btang cleaning"
  exit 1
fi

# TODO samtools flagstat logged on each bam?

# setup inputs
# mkdir 00_inputs
<%=
  #link_fastq_inputs()
%>

# ln -s ${FASTQ1} 00_inputs/A_1.fastq
# ln -s ${FASTQ2} 00_inputs/A_2.fastq
# ln -s ${FASTQ3} 00_inputs/B_1.fastq
# ln -s ${FASTQ4} 00_inputs/B_2.fastq

mkdir 01_bwa_aln_sai
# prep all reads for alignment
qsub -o logs -sync y -t 1-<%= total_number_input_sequence_files() %> -b y -V -j y -cwd -q all.q -N <%= sample_name %>_bwa_aln bwa_aln_qsub_tasked.rb 01_bwa_aln_sai /Volumes/hts_core/Shared/homo_sapiens_36.1/chr_fixed/bwa_indexed/hg18.fa <%= ordered_fastq_inputs() %>

if [ "$?" -ne "0" ]; then
  echo -e "Failure with bwa sai"
  exit 1
fi


mkdir 02_bwa_alignment
# align two lanes
<%=
  bwa_aligment_command(sample_name,data)
%>

if [ "$?" -ne "0" ]; then
  echo -e "Failure with bwa alignment"
  exit 1
fi

mkdir 03_first_bam
# make some bams
qsub -o logs -sync y -t 1-<%= total_number_input_sequenced_lanes() %> -b y -V -j y -cwd -q all.q -N <%= sample_name %>_make_bam make_bam_qsub_tasked.rb 03_first_bam ${GATK_REF} <%= input_sam_bam_files("02_bwa_alignment","sam") %>

if [ "$?" -ne "0" ]; then
  echo -e "Failure making first bams"
  exit 1
fi

mkdir 04_merged_bam
# merge the two lanes into one bam
qsub -o logs -sync y -b y -V -j y -cwd -q all.q -N <%= sample_name %>_merge_bams picard MergeSamFiles <%= input_sam_bam_files("INPUT=./03_first_bam","bam") %> OUTPUT=./04_merged_bam/cleaned.bam USE_THREADING=True

if [ "$?" -ne "0" ]; then
  echo -e "Failure merging bams"
  exit 1
fi

qsub -o logs -sync y -b y -V -j y -cwd -q all.q -N <%= sample_name %>_sort_merged samtools sort ./04_merged_bam/cleaned.bam 04_merged_bam/cleaned-sorted

if [ "$?" -ne "0" ]; then
  echo -e "Failure sorting bams"
  exit 1
fi

rm ./04_merged_bam/cleaned.bam && mv ./04_merged_bam/cleaned-sorted.bam ./04_merged_bam/cleaned.bam

qsub -o logs -sync y -b y -V -j y -cwd -q all.q -N <%= sample_name %>_index_merged samtools index ./04_merged_bam/cleaned.bam ./04_merged_bam/cleaned.bai

if [ "$?" -ne "0" ]; then
  echo -e "Failure indexing bams"
  exit 1
fi

mkdir 05_dup_marked
# mark duplicates with picard
qsub -o logs -sync y -b y -V -j y -cwd -q all.q -N <%= sample_name %>_mark_dups picard MarkDuplicates INPUT=./04_merged_bam/cleaned.bam OUTPUT=./05_dup_marked/cleaned.bam METRICS_FILE=./05_dup_marked/mark_dups_metrics.txt

if [ "$?" -ne "0" ]; then
  echo -e "Failure with marking the duplicates"
  exit 1
fi

qsub -o logs -sync y -b y -V -j y -cwd -q all.q -N <%= sample_name %>_index_merged_dup samtools index ./05_dup_marked/cleaned.bam ./05_dup_marked/cleaned.bam.bai

if [ "$?" -ne "0" ]; then                                                                                                                                                                                            
 echo -e "Failure indexing duplicate marked bam"                                                                                                                                                                       
 exit 1                                                                                                                                                                                                              
fi 

mkdir 06_intervals
# calculate intervals for realignment
qsub -o logs -sync y -b y -V -j y -cwd -q all.q -N <%= sample_name %>_intervals gatk -T RealignerTargetCreator -R ${GATK_REF} -I ./05_dup_marked/cleaned.bam -o ./06_intervals/cleaned.intervals

if [ "$?" -ne "0" ]; then
 echo -e "Failure with target realigment creation"
 exit 1
fi

mkdir 07_realigned_bam
# realign
qsub -o logs -sync y -b y -V -j y -cwd -q all.q -N <%= sample_name %>_realign gatk -T IndelRealigner -R ${GATK_REF} -I ./05_dup_marked/cleaned.bam --targetIntervals ./06_intervals/cleaned.intervals -o ./07_realigned_bam/cleaned.bam -maxInRam 1000000

if [ "$?" -ne "0" ]; then
 echo -e "Failure with indel realigmnent"
 exit 1
fi

qsub -o logs -sync y -b y -V -j y -cwd -q all.q -N <%= sample_name %>_fixmates picard FixMateInformation INPUT=./07_realigned_bam/cleaned.bam SO=coordinate VALIDATION_STRINGENCY=SILENT

if [ "$?" -ne "0" ]; then
  echo -e "Failure fixing mate info"
  exit 1
fi

qsub -o logs -sync y -b y -V -j y -cwd -q all.q -N <%= sample_name %>_index_fixed samtools index ./07_realigned_bam/cleaned.bam ./07_realigned_bam/cleaned.bai

if [ "$?" -ne "0" ]; then
  echo -e "Failure indexing"
  exit 1
fi


mkdir 08_uncalibated_covariates
# recalibration 
qsub -o logs -sync y -b y -V -j y -cwd -q all.q -N <%= sample_name %>_uncalibrated_covariates gatk -T CountCovariates -R ${GATK_REF} -D ${GATK_DBSNP} -I ./07_realigned_bam/cleaned.bam -cov ReadGroupCovariate -cov QualityScoreCovariate -cov CycleCovariate -cov DinucCovariate -recalFile ./08_uncalibated_covariates/recal_data.csv -nt 8

if [ "$?" -ne "0" ]; then
 echo -e "Failure counting covariates"
 exit 1
fi

mkdir 09_original_covariate_analysis
qsub -o logs -b y -V -j y -cwd -q all.q -N <%= sample_name %>_analyze_covariates java -Xmx4g -jar ${GATK_BASE}/resources/AnalyzeCovariates.jar -resources ${GATK_BASE}/resources -recalFile ./08_uncalibated_covariates/recal_data.csv -outputDir ./09_original_covariate_analysis

mkdir 10_recalibrated_bam
qsub -o logs -sync y -b y -V -j y -cwd -q all.q -N <%= sample_name %>_recalibrate gatk -T TableRecalibration -R ${GATK_REF} -I ./07_realigned_bam/cleaned.bam -recalFile ./08_uncalibated_covariates/recal_data.csv -o ./10_recalibrated_bam/recalibrated.bam <%= default_rg(sample_name,data) %>

if [ "$?" -ne "0" ]; then
 echo -e "Failure reclibrating bam"
 exit 1
fi

qsub -o logs -sync y -b y -V -j y -cwd -q all.q -N <%= sample_name %>_sort_recalibrated samtools sort ./10_recalibrated_bam/recalibrated.bam ./10_recalibrated_bam/recalibrated-sorted

if [ "$?" -ne "0" ]; then
  echo -e "Failure sorting recalibrated"
  exit 1
fi

rm ./10_recalibrated_bam/recalibrated.bam && mv ./10_recalibrated_bam/recalibrated-sorted.bam ./10_recalibrated_bam/recalibrated.bam

qsub -o logs -sync y -b y -V -j y -cwd -q all.q -N <%= sample_name %>_recalibated_realigned samtools index ./10_recalibrated_bam/recalibrated.bam ./10_recalibrated_bam/recalibrated.bam.bai


mkdir 11_calibated_covariates
qsub -o logs -sync y -b y -V -j y -cwd -q all.q -N <%= sample_name %>_calibrated_covariates gatk -T CountCovariates -R ${GATK_REF} -D ${GATK_DBSNP} -I ./10_recalibrated_bam/recalibrated.bam -cov ReadGroupCovariate -cov QualityScoreCovariate -cov CycleCovariate -cov DinucCovariate -recalFile ./11_calibated_covariates/recal_data.csv -nt 8

if [ "$?" -ne "0" ]; then
 echo -e "Failure counting calibrated covariates"
 exit 1
fi

mkdir 12_recalibrated_covariate_analysis
qsub -o logs -b y -V -j y -cwd -q all.q -N <%= sample_name %>_analyze_calibrated_covariates java -Xmx4g -jar ${GATK_BASE}/resources/AnalyzeCovariates.jar -resources ${GATK_BASE}/resources -recalFile ./11_calibated_covariates/recal_data.csv -outputDir ./12_recalibrated_covariate_analysis

mkdir 13_final_bam
# resort & index that bam
qsub -o logs -sync y -b y -V -j y -cwd -q all.q -N <%= sample_name %>_final_bam_sort samtools sort ./10_recalibrated_bam/recalibrated.bam ./13_final_bam/<%= sample_name %>
qsub -o logs -sync y -b y -V -j y -cwd -q all.q -N <%= sample_name %>_final_bam_index samtools index ./13_final_bam/<%= sample_name %>.bam ./13_final_bam/<%= sample_name %>.bam.bai

# fastqc
mkdir qc
qsub -o logs -b y -V -j y -cwd -q all.q -N <%= sample_name %>_qc fastqc -o qc <%= fastq_shell_vars() %> ./13_final_bam/<%= sample_name %>.bam

# call indels & snps
qsub -o logs -b y -V -j y -cwd -q all.q -N <%= sample_name %>_indels gatk -T IndelGenotyperV2 -R ${GATK_REF} -I ./13_final_bam/<%= sample_name %>.bam -o <%= sample_name %>_indels.vcf
qsub -o logs -sync y -b y -V -j y -cwd -q all.q -N <%= sample_name %>_snps gatk -T UnifiedGenotyper -R ${GATK_REF} -I ./13_final_bam/<%= sample_name %>.bam -o <%= sample_name %>_snps.vcf -stand_call_conf 30.0 -stand_emit_conf 10.0

if [ "$?" -ne "0" ]; then
 echo -e Failure
 exit 1
fi

rm -rf 00_inputs \
01_bwa_aln_sai \
02_bwa_alignment \
03_first_bam \
04_merged_bam \
05_dup_marked \
06_intervals \
07_realigned_bam \
08_uncalibated_covariates \
10_recalibrated_bam \
11_calibated_covariates

# gzip something
# TODO ${FASTQ1} ${FASTQ2} ${FASTQ3} ${FASTQ4} & rejects
<%=
  gzip_original_fastq(sample_name)
%>

touch finished.txt
EOF

# check them all first
samples.each do |s|
  raise "Can't find sample #{s}" unless config[s]
end

samples.each do |sample_name|
  data = config[sample_name]
  output_dir = File.join(output_base,sample_name)

  raise "Failed to make dir: #{output_dir} for #{sample_name}" unless Dir.mkdir(output_dir)

  script_file = File.join(output_dir,"analyze.sh")
  File.open(script_file,"w") do |f|
    f.puts ERB.new(script_template).result(binding)
  end
  
  return_dir = Dir.pwd
  raise "Failed to change to dir: #{output_dir} for #{sample_name}" unless Dir.chdir(output_dir)

  Dir.mkdir("logs")
  
  cmd = "qsub -o logs -sync y -b y -V -j y -cwd -q all.q -m e -N #{sample_name}_full ./analyze.sh"
  puts cmd
  system cmd
  
  Dir.chdir(return_dir)
end
