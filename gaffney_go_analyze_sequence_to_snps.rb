#!/usr/bin/env ruby1.9

require 'yaml'
require 'erb'

output_base = ARGV.shift
config = YAML::load(File.open(ARGV.shift))

samples = ARGV.clone

script_template = ERB.new <<-EOF
#!/bin/bash

GATK_REF=/Volumes/hts_core/Shared/homo_sapiens_36.1/hg18_with_mt.fasta
GATK_DBSNP=/Volumes/hts_core/Shared/dbsnp/dbsnp_129_hg18.rod

SAMPLE="lgs302805"

FASTQ1=`pwd`"/lgs302805_utsw_48_gagt_3/lgs302805_cleaned_utsw_48_gagt_3_1.fastq"
FASTQ2=`pwd`"/lgs302805_utsw_48_gagt_3/lgs302805_cleaned_utsw_48_gagt_3_2.fastq"
FASTQ3=`pwd`"/lgs302805_utsw_48_gagt_4/lgs302805_cleaned_utsw_48_gagt_4_1.fastq"
FASTQ4=`pwd`"/lgs302805_utsw_48_gagt_4/lgs302805_cleaned_utsw_48_gagt_4_2.fastq"

module load bwa/0.5.9rc1
module load samtools/0.1.12
module load picard/1.36
module load gatk/1.0.4705
module load fastqc/0.7.2
module load btangs/1.2.0

GATK_BIN=`which gatk`
GATK_BASE=`dirname ${GATK_BIN}`"/.."

mkdir logs

# clean
# TODO
clean_sample.rb -r ${RUN} -l ${LANE} -s ${SAMPLE} -b . [--single] INPUT_SEQUENCE

# TODO samtools flagstat logged on each bam?

# setup inputs
mkdir 00_inputs
ln -s ${FASTQ1} 00_inputs/A_1.fastq
ln -s ${FASTQ2} 00_inputs/A_2.fastq
ln -s ${FASTQ3} 00_inputs/B_1.fastq
ln -s ${FASTQ4} 00_inputs/B_2.fastq

mkdir 01_bwa_aln_sai
# prep all reads for alignment
qsub -o logs -sync y -t 1-4 -b y -V -j y -cwd -q all.q -N ${SAMPLE}_bwa_aln ~/tmp/bwa_aln_qsub_tasked.rb 01_bwa_aln_sai /Volumes/hts_core/Shared/homo_sapiens_36.1/bwa_indexed/hg18.fa 00_inputs/A_1.fastq 00_inputs/A_2.fastq 00_inputs/B_1.fastq 00_inputs/B_2.fastq

if [ "$?" -ne "0" ]; then
  echo -e Failure
  exit 1
fi


mkdir 02_bwa_sampe
# align two lanes
qsub -o logs -sync y -t 1-2 -b y -V -j y -cwd -q all.q -N ${SAMPLE}_bwa_sampe ~/tmp/bwa_sampe_qsub_tasked.rb 02_bwa_sampe /Volumes/hts_core/Shared/homo_sapiens_36.1/bwa_indexed/hg18.fa '"@RG\\tID:lgs300889_s_5\\tSM:lgs300889\\tPL:Illumina\\tPU:s_5"' 01_bwa_aln_sai/A_1.sai 01_bwa_aln_sai/A_2.sai 00_inputs/A_1.fastq 00_inputs/A_2.fastq '"@RG\\tID:lgs300889_s_6\\tSM:lgs300889\\tPL:Illumina\\tPU:s_6"' 01_bwa_aln_sai/B_1.sai 01_bwa_aln_sai/B_2.sai 00_inputs/B_1.fastq 00_inputs/B_2.fastq

if [ "$?" -ne "0" ]; then
  echo -e Failure
  exit 1
fi

mkdir 03_first_bam
# make some bams
qsub -o logs -sync y -t 1-2 -b y -V -j y -cwd -q all.q -N ${SAMPLE}_make_bam ~/tmp/make_bam_qsub_tasked.rb 03_first_bam ${GATK_REF} 02_bwa_sampe/A.sam 02_bwa_sampe/B.sam

if [ "$?" -ne "0" ]; then
  echo -e Failure
  exit 1
fi

mkdir 04_merged_bam
# merge the two lanes into one bam
qsub -o logs -sync y -b y -V -j y -cwd -q all.q -N ${SAMPLE}_merge_bams picard MergeSamFiles INPUT=./03_first_bam/A.bam INPUT=./03_first_bam/B.bam OUTPUT=./04_merged_bam/cleaned.bam USE_THREADING=True

if [ "$?" -ne "0" ]; then
  echo -e Failure
  exit 1
fi

qsub -o logs -sync y -b y -V -j y -cwd -q all.q -N ${SAMPLE}_sort_merged samtools sort ./04_merged_bam/cleaned.bam 04_merged_bam/cleaned-sorted

if [ "$?" -ne "0" ]; then
  echo -e Failure
  exit 1
fi

rm ./04_merged_bam/cleaned.bam && mv ./04_merged_bam/cleaned-sorted.bam ./04_merged_bam/cleaned.bam

qsub -o logs -sync y -b y -V -j y -cwd -q all.q -N ${SAMPLE}_index_merged samtools index ./04_merged_bam/cleaned.bam ./04_merged_bam/cleaned.bai

if [ "$?" -ne "0" ]; then
  echo -e Failure
  exit 1
fi

mkdir 05_dup_marked
# mark duplicates with picard
qsub -o logs -sync y -b y -V -j y -cwd -q all.q -N ${SAMPLE}_mark_dups picard MarkDuplicates INPUT=./04_merged_bam/cleaned.bam OUTPUT=./05_dup_marked/cleaned.bam METRICS_FILE=./05_dup_marked/mark_dups_metrics.txt

if [ "$?" -ne "0" ]; then
  echo -e Failure
  exit 1
fi

qsub -o logs -sync y -b y -V -j y -cwd -q all.q -N ${SAMPLE}_index_merged_dup samtools index ./05_dup_marked/cleaned.bam ./05_dup_marked/cleaned.bam.bai

if [ "$?" -ne "0" ]; then                                                                                                                                                                                            
 echo -e Failure                                                                                                                                                                                                     
 exit 1                                                                                                                                                                                                              
fi 

mkdir 06_intervals
# calculate intervals for realignment
qsub -o logs -sync y -b y -V -j y -cwd -q all.q -N ${SAMPLE}_intervals gatk -T RealignerTargetCreator -R ${GATK_REF} -I ./05_dup_marked/cleaned.bam -o ./06_intervals/cleaned.intervals

if [ "$?" -ne "0" ]; then
 echo -e Failure
 exit 1
fi

mkdir 07_realigned_bam
# realign
qsub -o logs -sync y -b y -V -j y -cwd -q all.q -N ${SAMPLE}_realign gatk -T IndelRealigner -R ${GATK_REF} -I ./05_dup_marked/cleaned.bam --targetIntervals ./06_intervals/cleaned.intervals -o ./07_realigned_bam/cleaned.bam -maxInRam 1000000

if [ "$?" -ne "0" ]; then
 echo -e Failure
 exit 1
fi

qsub -o logs -sync y -b y -V -j y -cwd -q all.q -N ${SAMPLE}_fixmates picard FixMateInformation INPUT=./07_realigned_bam/cleaned.bam SO=coordinate VALIDATION_STRINGENCY=SILENT

if [ "$?" -ne "0" ]; then
  echo -e Failure
  exit 1
fi

qsub -o logs -sync y -b y -V -j y -cwd -q all.q -N ${SAMPLE}_index_fixed samtools index ./07_realigned_bam/cleaned.bam ./07_realigned_bam/cleaned.bai

if [ "$?" -ne "0" ]; then
  echo -e Failure
  exit 1
fi


mkdir 08_uncalibated_covariates
# recalibration 
qsub -o logs -sync y -b y -V -j y -cwd -q all.q -N ${SAMPLE}_uncalibrated_covariates gatk -T CountCovariates -R ${GATK_REF} -D ${GATK_DBSNP} -I ./07_realigned_bam/cleaned.bam -cov ReadGroupCovariate -cov QualityScoreCovariate -cov CycleCovariate -cov DinucCovariate -recalFile ./08_uncalibated_covariates/recal_data.csv -nt 8

if [ "$?" -ne "0" ]; then
 echo -e Failure
 exit 1
fi

mkdir 09_original_covariate_analysis
qsub -o logs -b y -V -j y -cwd -q all.q -N ${SAMPLE}_analyze_covariates java -Xmx4g -jar ${GATK_BASE}/resources/AnalyzeCovariates.jar -resources ${GATK_BASE}/resources -recalFile ./08_uncalibated_covariates/recal_data.csv -outputDir ./09_original_covariate_analysis

if [ "$?" -ne "0" ]; then
 echo -e Failure
 exit 1
fi

mkdir 10_recalibrated_bam
qsub -o logs -sync y -b y -V -j y -cwd -q all.q -N ${SAMPLE}_recalibrate gatk -T TableRecalibration -R ${GATK_REF} -I ./07_realigned_bam/cleaned.bam -recalFile ./08_uncalibated_covariates/recal_data.csv -o ./10_recalibrated_bam/recalibrated.bam

if [ "$?" -ne "0" ]; then
 echo -e Failure
 exit 1
fi

qsub -o logs -sync y -b y -V -j y -cwd -q all.q -N ${SAMPLE}_sort_recalibrated samtools sort ./10_recalibrated_bam/recalibrated.bam ./10_recalibrated_bam/recalibrated-sorted

if [ "$?" -ne "0" ]; then
  echo -e Failure
  exit 1
fi

rm ./10_recalibrated_bam/recalibrated.bam && mv ./10_recalibrated_bam/recalibrated-sorted.bam ./10_recalibrated_bam/recalibrated.bam

qsub -o logs -sync y -b y -V -j y -cwd -q all.q -N ${SAMPLE}_recalibated_realigned samtools index ./10_recalibrated_bam/recalibrated.bam ./10_recalibrated_bam/recalibrated.bam.bai


mkdir 11_calibated_covariates
qsub -o logs -sync y -b y -V -j y -cwd -q all.q -N ${SAMPLE}_calibrated_covariates gatk -T CountCovariates -R ${GATK_REF} -D ${GATK_DBSNP} -I ./10_recalibrated_bam/recalibrated.bam -cov ReadGroupCovariate -cov QualityScoreCovariate -cov CycleCovariate -cov DinucCovariate -recalFile ./11_calibated_covariates/recal_data.csv -nt 8

if [ "$?" -ne "0" ]; then
 echo -e Failure
 exit 1
fi

mkdir 12_recalibrated_covariate_analysis
qsub -o logs -b y -V -j y -cwd -q all.q -N ${SAMPLE}_analyze_calibrated_covariates java -Xmx4g -jar ${GATK_BASE}/resources/AnalyzeCovariates.jar -resources ${GATK_BASE}/resources -recalFile ./11_calibated_covariates/recal_data.csv -outputDir ./12_recalibrated_covariate_analysis

if [ "$?" -ne "0" ]; then
 echo -e Failure
 exit 1
fi

mkdir 13_final_bam
# resort & index that bam
qsub -o logs -sync y -b y -V -j y -cwd -q all.q -N ${SAMPLE}_final_bam_sort samtools sort ./10_recalibrated_bam/recalibrated.bam ./13_final_bam/${SAMPLE}
qsub -o logs -sync y -b y -V -j y -cwd -q all.q -N ${SAMPLE}_final_bam_index samtools index ./13_final_bam/${SAMPLE}.bam ./13_final_bam/${SAMPLE}.bam.bai

# fastqc
mkdir qc
qsub -o logs -b y -V -j y -cwd -q all.q -N ${SAMPLE}_qc fastqc -o qc ${FASTQ1} ${FASTQ2} ${FASTQ3} ${FASTQ4} ./13_final_bam/${SAMPLE}.bam

# call indels & snps
qsub -o logs -b y -V -j y -cwd -q all.q -N ${SAMPLE}_indels gatk -T IndelGenotyperV2 -R ${GATK_REF} -I ./13_final_bam/${SAMPLE}.bam -o ${SAMPLE}_indels.vcf
qsub -o logs -sync y -b y -V -j y -cwd -q all.q -N ${SAMPLE}_snps gatk -T UnifiedGenotyper -R ${GATK_REF} -I ./13_final_bam/${SAMPLE}.bam -o ${SAMPLE}_snps.vcf -stand_call_conf 30.0 -stand_emit_conf 10.0 -nt 8

if [ "$?" -ne "0" ]; then
 echo -e Failure
 exit 1
fi

rm -rf 00_inputs \
01_bwa_aln_sai \
02_bwa_sampe \
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
qsub -o logs -b y -V -j y -cwd -q all.q -N ${SAMPLE}_gzip_1 gzip --fast ${FASTQ1}
qsub -o logs -b y -V -j y -cwd -q all.q -N ${SAMPLE}_gzip_1 gzip --fast ${FASTQ2}
qsub -o logs -b y -V -j y -cwd -q all.q -N ${SAMPLE}_gzip_1 gzip --fast ${FASTQ3}
qsub -o logs -b y -V -j y -cwd -q all.q -N ${SAMPLE}_gzip_1 gzip --fast ${FASTQ4}
EOF

# check them all first
samples.each do |s|
  raise "Can't find sample #{s}" unless config[s]
end

samples.each do |s|
  data = config[s]
  output_dir = File.join(output_base,s)

  raise "Failed to make dir: #{output_dir} for #{s}" unless Dir.mkdir(output_dir)

  script_file = File.join(output_dir,"analyze.sh")
  File.open(script_file,"w") do |f|
    f.puts script_template.result(binding)
  end
  
  return_dir = Dir.pwd
  raise "Failed to change to dir: #{output_dir} for #{s}" unless Dir.chdir(output_dir)

  cmd = "qsub -o logs -sync y -b y -V -j y -cwd -q all.q -m e -N ${SAMPLE}_full #{script_file}"
  puts cmd
  system cmd
  
  Dir.chdir(return_dir)
end