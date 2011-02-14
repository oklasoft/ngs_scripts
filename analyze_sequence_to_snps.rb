#!/usr/bin/env ruby1.9
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
#  Copyright (c) 2011 Stuart Glenn, Oklahoma Medical Research Foundation. (OMRF)
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


class AnalysisTemplate
  
  # Create instance
  # *+default_config+ - Default options
  # *+sample_name+ - Name of this sample
  # *+data+ - Hash of options for this sample, overrides anything in default_config
  def initialize(default_config,sample_name,data)
    @default_config = default_config
    @sample_name = sample_name
    @data = data
  end #initialize(default_config,sample_name,data)
  
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
        path = "`pwd`\"/#{prefix}/#{base_file}.fastq\""
        fastqs << "#{shell_var}=#{path}"
        @fastq_shell_vars[shell_var] = {:path  => path, :paired => pair_part, :letter => letter, :base_file => base_file, :prefix => prefix}
        @fastq_shell_vars_by_lane[-1] << shell_var
        @input_sam_files << {:index => s_i, :b_index => pair_part}
      end
    end
    fastqs.join("\n")
  end

  def fastq_shell_vars()
    @fastq_shell_vars.keys.map{|v| "${#{v}}"}.join(" ")
  end


  def ordered_sam_inputs()
    @input_sam_files.map {|s| "#{s[:b_index]} ./00_inputs/#{s[:index]}.bam"}.join(" ")
  end

  def cleanup_cleaned_fastq_files(sample_name)
    #qsub -o logs -b y -V -j y -cwd -q all.q -N <%= sample_name %>_gzip_1 gzip --fast ${FASTQ1}
    cmds = []
    @fastq_shell_vars_by_lane.flatten.each_with_index do |input,i|
      cmds << "rm -f ${#{input}}"
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

  def bwa_reference_for_data(data)
    data.first[:bwa_ref] || @default_config[:bwa_ref]
  end

  def bwa_alignment_command(sample_name,data)
    cmd = "qsub -o logs -sync y -t 1-#{total_number_input_sequenced_lanes()} -b y -V -j y -cwd -q all.q -N #{sample_name}_bwa_alignment bwa_sampese_qsub_tasked.rb 02_bwa_alignment #{bwa_reference_for_data(data)}"
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
        # cmd += " 01_bwa_aln_sai/#{@fastq_shell_vars[v][:base_file]}.sai"
        cmd += " 01_bwa_aln_sai/#{index}-#{@fastq_shell_vars[v][:paired]}.sai"
      end
      # fastq file(s)
      lane_shell_vars.each do |v|
        # cmd += " ${#{v}}"
        cmd += " 00_inputs/#{index}.bam"
      end
    end
    return cmd  
  end

  def create_original_sam_inputs(sample_name,data)
    cmd = "qsub -o logs -sync y -t 1-#{total_number_input_sequenced_lanes()} -b y -V -j y -cwd -q all.q -N #{sample_name}_convert_to_sam convert_fastq_to_sam_qsub_tasked.rb 00_inputs Illumina #{sample_name}"
    @fastq_shell_vars_by_lane.each_with_index do |lane_shell_vars,index|
      if data[index][:is_paired]
        cmd += " paired"
      else
        cmd += " single"
      end

      cmd += " #{sample_name}_#{data[index][:run]}_s_#{data[index][:lane]} #{data[index][:lane]}"

      # fastq file(s)
      lane_shell_vars.each do |v|
        cmd += " ${#{v}}"
      end
    end
    return cmd  
  end

  def fix_sam_read_group(sample_name,data)
    cmd = "qsub -o logs -sync y -t 1-#{total_number_input_sequenced_lanes()} -b y -V -j y -cwd -q all.q -N #{sample_name}_fix_sam_rg fix_sam_rg_qsub_tasked.rb 02_bwa_alignment Illumina #{sample_name}"
    @fastq_shell_vars_by_lane.each_with_index do |lane_shell_vars,index|
      # rg
      cmd += " #{sample_name}_#{data[index][:run]}_s_#{data[index][:lane]}"
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

  def covariate_or_final(sample_name,data)
    if data.first[:snp_rod] || @default_config[:snp_rod]
      covariate_recalibration(sample_name,data)
    else
      skip_covariate_recalibration(sample_name,data)
    end
  end

  def skip_covariate_recalibration(sample_name,data)
    <<-EOF
    mkdir 13_final_bam
    mv ./07_realigned_bam/cleaned.bam ./13_final_bam/#{sample_name}.bam
    mv ./07_realigned_bam/cleaned.bai ./13_final_bam/#{sample_name}.bai
    EOF
  end

  def covariate_recalibration(sample_name,data)
    ERB.new(<<-EOF
    # HERE

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
  EOF
    ).result(binding)
  end

  def script_template()
    return DATA.read
  end

  def to_s
    ERB.new(script_template()).result(binding)
  end
end


class AnalysisTemplaterApp
  VERSION       = "0.0.9-pre01"
  REVISION_DATE = "2011-01-10"
  AUTHOR        = "Stuart Glenn <Stuart-Glenn@omrf.org>"
  COPYRIGHT     = "Copyright (c) 2011 Oklahoma Medical Research Foundation"
  
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
      return 0
    end
    return run_real()
  end #run

  private
  
  # Do the work of running this app
  def run_real
    @default_config = (@config["DEFAULT"] || [{:run => nil, :lane => nil, :bwa_ref => nil, :gatk_ref => nil, :snp_rod => nil}]).first
    
    statii = @options.samples.map do |sample_name|
      process_sample(sample_name)
    end
    return statii.max
  end #run
  
  # Make the dir, the analyze script and launch said script for sample_name
  def process_sample(sample_name)
    data = @config[sample_name]
    output_dir = File.join(@options.output_base,sample_name)

    unless Dir.mkdir(output_dir)
      @stderr.puts "Failed to make dir: #{output_dir} for #{sample_name}" 
      return 1
    end

    script_file = File.join(output_dir,"analyze.sh")
    File.open(script_file,"w") do |f|
      f.puts AnalysisTemplate.new(@default_config,sample_name,data)
    end

    return_dir = Dir.pwd
    unless Dir.chdir(output_dir)
      @stderr.puts "Failed to change to dir: #{output_dir} for #{sample_name}" 
      return 1
    end

    Dir.mkdir("logs")

    # We sleep a random amount to avoid overloading SGE with a billion jobs right away
    sleep(rand(@options.delay))
    cmd = "qsub -o logs -sync y -b y -V -j y -cwd -q all.q -m e -N #{sample_name}_full ./analyze.sh"
    @stdout.puts(cmd)
    system(cmd)
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

      o.on("-d","--delay", "=REQUIRED") do |amount|
        @options.delay = amount.to_i
      end

      o.on("-c","--config", "=REQUIRED") do |conf_file|
        @options.config_file = conf_file
      end

      o.on("-o","--output", "=REQUIRED") do |output_destination|
        @options.output_base = output_destination
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
      :delay => 30
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
#!/bin/bash

# be sure to start with a fresh & known enviroment (hopefully)
module unload bwa
module load bwa/0.5.9
module unload samtools
module load samtools/0.1.12
module unload picard
module load picard/1.36
module unload gatk
module load gatk/1.0.4905
module unload fastqc
module load fastqc/0.7.2
module unload tabix
module load tabix/0.2.3
module unload btangs
module load btangs/1.2.0

# It is easier to use variables for thse paths
GATK_REF=<%= @data.first[:gatk_ref] || @default_config[:gatk_ref] %>
GATK_DBSNP=<%= @data.first[:snp_rod] || @default_config[:snp_rod] || 'nil' %>
GATK_BIN=`which gatk`
GATK_BASE=`dirname ${GATK_BIN}`"/.."

SAMPLE="<%= @sample_name %>"

<%=
  fastq_file_list(@sample_name,@data)
%>

# initial PCR clean by 'bins' via btangs
<%=
  clean_commands(@sample_name,@data)
%>

if [ "$?" -ne "0" ]; then
  echo -e "Failure with btang cleaning"
  exit 1
fi

# TODO samtools flagstat logged on each bam?

# setup input sams, will get illumina scores to standard sanger
mkdir 00_inputs
<%=
  create_original_sam_inputs(@sample_name,@data)
%>

if [ "$?" -ne "0" ]; then
  echo -e "Failure converting fastq to sam"
  exit 1
fi

# Prep all those reads for alignment with bwa aln 
mkdir 01_bwa_aln_sai
qsub -o logs -sync y -t 1-<%= total_number_input_sequence_files() %> -b y -V -j y -cwd -q all.q -N <%= @sample_name %>_bwa_aln bwa_aln_qsub_tasked.rb 01_bwa_aln_sai <%= bwa_reference_for_data(@data) %> <%= ordered_sam_inputs() %>

if [ "$?" -ne "0" ]; then
  echo -e "Failure with bwa sai"
  exit 1
fi

# Now take those & actually generate the alignment SAM output, paired or single
mkdir 02_bwa_alignment
<%=
  bwa_alignment_command(@sample_name,@data)
%>

if [ "$?" -ne "0" ]; then
  echo -e "Failure with bwa alignment"
  exit 1
fi

# Since bwa right now doesn't add read group info correctly (it is not doing it for unmapped), manually add read group
<%=
  # fix_sam_read_group(@sample_name,@data)
%>

# if [ "$?" -ne "0" ]; then
#   echo -e "Failure with SAM RG fix"
#   exit 1
# fi

# Going forward let us play with BAM files
mkdir 03_first_bam
qsub -o logs -sync y -t 1-<%= total_number_input_sequenced_lanes() %> -b y -V -j y -cwd -q all.q -N <%= @sample_name %>_make_bam make_bam_qsub_tasked.rb 03_first_bam ${GATK_REF} <%= input_sam_bam_files("02_bwa_alignment","sam") %>

if [ "$?" -ne "0" ]; then
  echo -e "Failure making first bams"
  exit 1
fi

# Now we might have had many input sets, so let us merge those all into a single BAM using picard
# TODO be smarter if there was only a single input
mkdir 04_merged_bam
qsub -o logs -sync y -b y -V -j y -cwd -q all.q -N <%= @sample_name %>_merge_bams picard MergeSamFiles <%= input_sam_bam_files("INPUT=./03_first_bam","bam") %> OUTPUT=./04_merged_bam/cleaned.bam USE_THREADING=True VALIDATION_STRINGENCY=LENIENT

if [ "$?" -ne "0" ]; then
  echo -e "Failure merging bams"
  exit 1
fi

# Make sure it really is sorted & indexed, sometimes things like to complain, just don't let them
qsub -o logs -sync y -b y -V -j y -cwd -q all.q -N <%= @sample_name %>_sort_merged samtools sort ./04_merged_bam/cleaned.bam 04_merged_bam/cleaned-sorted

if [ "$?" -ne "0" ]; then
  echo -e "Failure sorting bams"
  exit 1
fi

rm ./04_merged_bam/cleaned.bam && mv ./04_merged_bam/cleaned-sorted.bam ./04_merged_bam/cleaned.bam

qsub -o logs -sync y -b y -V -j y -cwd -q all.q -N <%= @sample_name %>_index_merged samtools index ./04_merged_bam/cleaned.bam ./04_merged_bam/cleaned.bai

if [ "$?" -ne "0" ]; then
  echo -e "Failure indexing bams"
  exit 1
fi

# Mark duplicates with picard
mkdir 05_dup_marked
qsub -o logs -sync y -b y -V -j y -cwd -q all.q -N <%= @sample_name %>_mark_dups picard MarkDuplicates INPUT=./04_merged_bam/cleaned.bam OUTPUT=./05_dup_marked/cleaned.bam METRICS_FILE=./05_dup_marked/mark_dups_metrics.txt VALIDATION_STRINGENCY=LENIENT

if [ "$?" -ne "0" ]; then
  echo -e "Failure with marking the duplicates"
  exit 1
fi

qsub -o logs -sync y -b y -V -j y -cwd -q all.q -N <%= @sample_name %>_index_merged_dup samtools index ./05_dup_marked/cleaned.bam ./05_dup_marked/cleaned.bam.bai

if [ "$?" -ne "0" ]; then                                                                                                                                                                                            
 echo -e "Failure indexing duplicate marked bam"                                                                                                                                                                       
 exit 1                                                                                                                                                                                                              
fi 

# Calculate intervals for realignment
mkdir 06_intervals
qsub -o logs -sync y -b y -V -j y -cwd -q all.q -N <%= @sample_name %>_intervals gatk -T RealignerTargetCreator -R ${GATK_REF} -I ./05_dup_marked/cleaned.bam -o ./06_intervals/cleaned.intervals

if [ "$?" -ne "0" ]; then
 echo -e "Failure with target realigment creation"
 exit 1
fi

# Now realign & fix any mate info
mkdir 07_realigned_bam
qsub -o logs -sync y -b y -V -j y -cwd -q all.q -N <%= @sample_name %>_realign gatk -T IndelRealigner -R ${GATK_REF} -I ./05_dup_marked/cleaned.bam --targetIntervals ./06_intervals/cleaned.intervals -o ./07_realigned_bam/cleaned.bam -maxInRam 1000000

if [ "$?" -ne "0" ]; then
 echo -e "Failure with indel realigmnent"
 exit 1
fi

qsub -o logs -sync y -b y -V -j y -cwd -q all.q -N <%= @sample_name %>_fixmates picard FixMateInformation INPUT=./07_realigned_bam/cleaned.bam SO=coordinate VALIDATION_STRINGENCY=SILENT

if [ "$?" -ne "0" ]; then
  echo -e "Failure fixing mate info"
  exit 1
fi

qsub -o logs -sync y -b y -V -j y -cwd -q all.q -N <%= @sample_name %>_index_fixed samtools index ./07_realigned_bam/cleaned.bam ./07_realigned_bam/cleaned.bai

if [ "$?" -ne "0" ]; then
  echo -e "Failure indexing"
  exit 1
fi

<%= covariate_or_final(@sample_name,@data) %>

# fastqc info
mkdir qc
qsub -o logs -b y -V -j y -cwd -q all.q -N <%= @sample_name %>_qc fastqc -o qc <%= fastq_shell_vars() %> ./13_final_bam/<%= @sample_name %>.bam

# Finally call individuals indels & snps
qsub -o logs -b y -V -j y -cwd -q all.q -N <%= @sample_name %>_indels gatk -T UnifiedGenotyper -glm DINDEL -R ${GATK_REF} -I ./13_final_bam/<%= @sample_name %>.bam -o <%= @sample_name %>_indels.vcf
qsub -o logs -sync y -b y -V -j y -cwd -q all.q -N <%= @sample_name %>_snps gatk -T UnifiedGenotyper -R ${GATK_REF} -I ./13_final_bam/<%= @sample_name %>.bam -o <%= @sample_name %>_snps.vcf -stand_call_conf 30.0 -stand_emit_conf 10.0

if [ "$?" -ne "0" ]; then
 echo -e Failure
 exit 1
fi

bgzip <%= @sample_name %>_snps.vcf
tabix <%= @sample_name %>_snps.vcf.gz

# Clean up after ourselves
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
11_calibated_covariates \
qc/*.zip

# gzip & clean up some of the cleaned input as we keep some it for now, but don't need it fullsized
<%=
  cleanup_cleaned_fastq_files(@sample_name)
%>

touch finished.txt