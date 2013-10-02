#!/usr/bin/env ruby1.9.3
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


class Template

  # Create instance
  # *+default_config+ - Default options
  # *+sample_name+ - Name of this sample
  # *+data+ - Hash of options for this sample, overrides anything in default_config
  def initialize(default_config,sample_name,data)
    @default_config = default_config
    @sample_name = sample_name
    @data = data
  end #initialize(default_config,sample_name,data)

  def to_s
    ERB.new(script_template()).result(binding)
  end

  def bash_header()
    <<-EOS
#!/bin/bash

cd "$( dirname "${BASH_SOURCE[0]}" )"

# be sure to start with a fresh & known enviroment (hopefully)
source /etc/profile.d/*.sh
module unload bwa
module load bwa/0.7.5a
module unload samtools
module load samtools/0.1.19
module unload picard
module load picard/1.99
module unload gatk
module load gatk/2.7-2-g6bda569
module unload fastqc
module load fastqc/0.10.1
module unload tabix
module load tabix/0.2.6
module unload btangs
module load btangs/1.6.0
    EOS
  end
end

class UnalignedExtractTemplate < Template

  def bam_to_fastq_for_run(run,bam_index)
    cmd = "bam2fastq --no-aligned --unaligned -o unaligned_fastq/#{@sample_name}_unaligned_#{run[:run]}_#{run[:lane]}\#_sequence.txt 03_first_bam/#{bam_index}.bam"
    cmd = "qsub -o logs -b y -j y -cwd -V -q ngs.q -l h_vmem=8G -sync y -N a_#{@sample_name}_unaligned_#{bam_index} #{cmd}"
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
 echo -e "Failure extracting unalinged from bam <%= i %>"
 exit 1
fi
<% end %>
qsub -o logs -b y -V -j y -cwd -q ngs.q -sync y -N a_<%=@sample_name%>_unaligned_gzip gzip -7 unaligned_fastq/*_sequence.txt

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
        path = if @default_config[:opts][:skip_btangs]
          input
        else
          "`pwd`\"/#{prefix}/#{base_file}.fastq\""
        end
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


  def ordered_bam_inputs()
    @input_sam_files.map {|s| "#{s[:b_index]} ./00_inputs/#{s[:index]}.bam"}.join(" ")
  end

  def cleanup_cleaned_fastq_files(sample_name)
    return "echo noop" if @default_config[:opts][:skip_btangs]
    #qsub -o logs -b y -V -j y -cwd -q ngs.q -N a_<%= sample_name %>_gzip_1 gzip --fast ${FASTQ1}
    cmds = []
    @fastq_shell_vars_by_lane.flatten.each_with_index do |input,i|
      cmds << "rm -f ${#{input}}"
      cmds << "qsub -o logs -b y -V -j y -cwd -q ngs.q -N a_#{sample_name}_gzip_#{i}_rejects gzip -9 #{@fastq_shell_vars[input][:prefix]}/rejects.txt" if 0==i%2
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
    cmd = "qsub -pe threaded 10 -l virtual_free=800M,h_vmem=32G -o logs -sync y -t 1-#{total_number_input_sequenced_lanes()} -b y -V -j y -cwd -q ngs.q -N a_#{sample_name}_bwa_alignment bwa_mem_qsub_tasked.rb 03_sorted_bams #{bwa_reference_for_data(data)}"
    @fastq_shell_vars_by_lane.each_with_index do |lane_shell_vars,index|
      if data[index][:is_paired]
        cmd += " paired"
      else
        cmd += " single"
      end

      # rg
      cmd += " '\"@RG\\\\tID:#{sample_name}_#{data[index][:run]}_s_#{data[index][:lane]}\\\\tSM:#{sample_name}\\\\tPL:Illumina\\\\tPU:#{data[index][:lane]}\"'"

      # 01_bwa_aln_sai file(s)
      #lane_shell_vars.each do |v|
        ## cmd += " 01_bwa_aln_sai/#{@fastq_shell_vars[v][:base_file]}.sai"
        #cmd += " 01_bwa_aln_sai/#{index}-#{@fastq_shell_vars[v][:paired]}.sai"
      #end
      # fastq file(s)
      lane_shell_vars.each do |v|
        cmd += " ${#{v}}"
        #cmd += " 00_inputs/#{index}.bam"
      end
    end
    return cmd
  end

  def extract_unaligned(sample_name,data)

  end

  def create_original_sam_inputs(sample_name,data)
    cmd = "qsub -l virtual_free=10G,h_vmem=50G -o logs -sync y -t 1-#{total_number_input_sequenced_lanes()} -b y -V -j y -cwd -q ngs.q -N a_#{sample_name}_convert_to_sam convert_fastq_to_bam_qsub_tasked.rb 00_inputs Illumina #{sample_name}"
    @fastq_shell_vars_by_lane.each_with_index do |lane_shell_vars,index|
      if data[index][:is_paired]
        cmd += " paired"
      else
        cmd += " single"
      end

      cmd += " #{sample_name}_#{data[index][:run]}_s_#{data[index][:lane]} #{data[index][:lane]} #{data[index][:quality_type] || "Illumina"}"

      # fastq file(s)
      lane_shell_vars.each do |v|
        cmd += " ${#{v}}"
      end
    end
    return cmd
  end

  def fix_sam_read_group(sample_name,data)
    cmd = "qsub -o logs -sync y -t 1-#{total_number_input_sequenced_lanes()} -b y -V -j y -cwd -q ngs.q -N a_#{sample_name}_fix_sam_rg fix_sam_rg_qsub_tasked.rb 02_bwa_alignment Illumina #{sample_name}"
    @fastq_shell_vars_by_lane.each_with_index do |lane_shell_vars,index|
      # rg
      cmd += " #{sample_name}_#{data[index][:run]}_s_#{data[index][:lane]}"
    end
    return cmd
  end

  def default_rg(sample_name,data)
    return ""
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
    return if @default_config[:opts][:skip_btangs]
    #clean_sample.rb -r ${RUN} -l ${LANE} -s <%= sample_name %> -b . [--single] INPUT_SEQUENCE
    cleans = []
    data.each_with_index do |sequence,s_i|
      cmd = "clean_sample.rb -s #{sample_name} -r #{sequence[:run]} -l #{sequence[:lane]} -b ."
      if sequence[:trim_end]
        cmd += " --trim-end #{sequence[:trim_end]}"
      end
      unless sequence[:is_paired]
        cmd += " --single-end"
      end
      cmd += " #{sequence[:inputs].join(" ").gsub(/\\/,"\\\\\\")}"
      cleans << "qsub -pe threaded 2 -l hadoop=1,h_vmem=4G -o logs -sync y -b y -V -j y -cwd -q ngs.q -N a_#{sample_name}_clean_#{s_i+1} #{cmd}"
    end
    cleans.join("\n") + <<-EOF

    if [ "$?" -ne "0" ]; then
      echo -e "Failure with btang cleaning"
      exit 1
    fi 

EOF
  end

  # output a -D SNP rod file, if we have a snp_rod
  def opt_d_rod_path(data)
    if data.first[:snp_rod] || @default_config[:snp_rod]
      "-D ${GATK_DBSNP}"
    else
      ""
    end
  end

  

  def variant_call(sample_name,data)
    if @default_config[:opts][:skip_vcf]
      return ""
    else
      bam_dir = "13_final_bam"
      bam_dir = "14_reduced_bam" if @default_config[:opts][:reduce_reads]
      ERB.new(<<-EOF
      # Finally call individuals indels & snps

      qsub -pe threaded 2 -o logs -sync y -b y -V -j y -cwd -q ngs.q -N a_<%= sample_name %>_variants -l mem_free=4G,h_vmem=6G gatk -T UnifiedGenotyper -A AlleleBalance -l INFO -nct 3 -R ${GATK_REF} -glm BOTH -I ./<%= bam_dir %>/<%= sample_name %>.bam -o <%= sample_name %>_variants.vcf -stand_call_conf <%= unified_genotyper_strand_call_conf(data) %> -stand_emit_conf 10.0 <%= opt_d_rod_path(data) %>

      if [ "$?" -ne "0" ]; then
       echo -e Failure
       exit 1
      fi

      bgzip <%= sample_name %>_variants.vcf
      tabix <%= sample_name %>_variants.vcf.gz
      EOF
      ).result(binding)
    end
  end

  def reduce_reads(sample_name,data)
    chr_bams = %w/1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y MT/.map do |chr|
      "INPUT=14_reduced_bam/by_chr/#{sample_name}-#{chr}.bam"
    end
    if @default_config[:opts][:reduce_reads]
      ERB.new(<<-EOF
      # Make a reduce reads BAM for variant calling better/faster/stronger
      mkdir 14_reduced_bam
      export JAVA_MEM_OPTS="-Xmx12G"
      qsub -t 1-25 -o logs -sync y -b y -V -j y -cwd -q ngs.q -N a_<%= sample_name %>_reduce_reads -l mem_free=8G,h_vmem=15G read_reducer_qsub_tasked.rb ${GATK_REF} 14_reduced_bam <%= sample_name %> ./13_final_bam/<%= sample_name %>.bam

      if [ "$?" -ne "0" ]; then
       echo -e "Failure reducing reads"
       exit 1
      fi

      export JAVA_MEM_OPTS="-Xmx20G"
      qsub -o logs -sync y -b y -V -j y -cwd -q ngs.q -N a_<%= sample_name %>_join_reduce_reads -l mem_free=12G,h_vmem=48G picard MergeSamFiles TMP_DIR=./tmp OUTPUT=14_reduced_bam/<%= sample_name %>.bam USE_THREADING=True VALIDATION_STRINGENCY=LENIENT MAX_RECORDS_IN_RAM=3000000 COMPRESSION_LEVEL=7 CREATE_INDEX=True SORT_ORDER=coordinate <%= chr_bams.join(" ") %>

      if [ "$?" -ne "0" ]; then
       echo -e "Failure joining reduced reads"
       exit 1
      fi
      mv 14_reduced_bam/<%= sample_name %>.bai 14_reduced_bam/<%= sample_name %>.bam.bai
      rm -rf 14_reduced_bam/by_chr

      EOF
      ).result(binding)
    else
      return ""
    end
  end
  
  def known_indels_opts()
    if @default_config[:known_indels]
      return @default_config[:known_indels].map {|i| "-known #{i}"}.join(" ")
    else
      return ""
    end
  end

  def indel_realign(sample_name,data)
    if @default_config[:opts][:skip_indel_realign]
      return <<-EOF
        ln -s 05_dup_marked 07_realigned_bam
      EOF
    else
      indel_realignment(sample_name,data)
    end
  end

  def indel_realignment(sample_name,data)
    compression = if data.first[:recalibration_known_sites] || @default_config[:recalibration_known_sites]
                    ""
                  else
                    "--bam_compression 7"
                  end
    ERB.new(<<-EOF
      # Calculate intervals for realignment
      mkdir 06_intervals
      qsub -pe threaded 10 -R y -o logs -sync y -b y -V -j y -cwd -q ngs.q -N a_<%= sample_name %>_intervals -l mem_free=1G,h_vmem=20G gatk -T RealignerTargetCreator -R ${GATK_REF} -I ./05_dup_marked/cleaned.bam -o ./06_intervals/cleaned.intervals -nt 10

      if [ "$?" -ne "0" ]; then
       echo -e "Failure with target realigment creation"
       exit 1
      fi

      # Now realign & fix any mate info
      mkdir 07_realigned_bam
      unset JAVA_MEM_OPTS
      qsub -o logs -sync y -b y -V -j y -cwd -q ngs.q -N a_<%= sample_name %>_realign -l mem_free=4G,h_vmem=6G gatk -T IndelRealigner <%= known_indels_opts() %> -R ${GATK_REF} -I ./05_dup_marked/cleaned.bam --targetIntervals ./06_intervals/cleaned.intervals -o ./07_realigned_bam/cleaned.bam --maxReadsInMemory 1000000 #{compression}

      if [ "$?" -ne "0" ]; then
       echo -e "Failure with indel realigmnent"
       exit 1
      fi
    EOF
    ).result(binding)
  end

  def mark_dupes_or_skip(sample_name,data)
    if @default_config[:opts][:skip_dupes]
      <<-EOF
qsub -l virtual_free=8G,h_vmem=56G -o logs -sync y -b y -V -j y -cwd -q ngs.q -N a_#{sample_name}_merge picard MergeSamFiles TMP_DIR=./tmp #{input_sam_bam_files("INPUT=03_sorted_bams","bam")} OUTPUT=./05_dup_marked/cleaned.bam VALIDATION_STRINGENCY=LENIENT MAX_RECORDS_IN_RAM=3000000 CREATE_INDEX=True USE_THREADING=True
if [ "$?" -ne "0" ]; then
  echo -e "Failure with merging the sams"
  exit 1
fi
      EOF
    else
      mark_dupes(sample_name,data)
    end
  end
  
  def mark_dupes(sample_name,data)
    <<-EOF
qsub -l virtual_free=8G,h_vmem=56G -o logs -sync y -b y -V -j y -cwd -q ngs.q -N a_#{sample_name}_merge_mark_dups picard MarkDuplicates TMP_DIR=./tmp #{input_sam_bam_files("INPUT=03_sorted_bams","bam")} OUTPUT=./05_dup_marked/cleaned.bam METRICS_FILE=./05_dup_marked/mark_dups_metrics.txt VALIDATION_STRINGENCY=LENIENT MAX_RECORDS_IN_RAM=3000000 CREATE_INDEX=True

if [ "$?" -ne "0" ]; then
  echo -e "Failure with marking the duplicates"
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
    mkdir 13_final_bam
    mv ./07_realigned_bam/cleaned.bam ./13_final_bam/#{sample_name}.bam
    mv ./07_realigned_bam/cleaned.bai ./13_final_bam/#{sample_name}.bai
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
    qsub -pe threaded 6 -R y -o logs -sync y -b y -V -j y -cwd -q ngs.q -N a_<%= sample_name %>_bqsr -l mem_free=4G,h_vmem=6G gatk -T BaseRecalibrator -R ${GATK_REF} <%= recalibration_known_sites() %> -I ./07_realigned_bam/cleaned.bam -o ./08_uncalibated_covariates/recal_data.grp -nct 8

    if [ "$?" -ne "0" ]; then
     echo -e "Failure counting covariates"
     exit 1
    fi

    mkdir 10_recalibrated_bam
    unset JAVA_MEM_OPTS
    qsub -pe threaded 6 -R y -o logs -sync y -b y -V -j y -cwd -q ngs.q -N a_<%= sample_name %>_recalibrate -l mem_free=4G,h_vmem=6G gatk -T PrintReads -R ${GATK_REF} -I ./07_realigned_bam/cleaned.bam -BQSR ./08_uncalibated_covariates/recal_data.grp -o ./10_recalibrated_bam/recalibrated.bam --bam_compression 7 -nct 8

    if [ "$?" -ne "0" ]; then
     echo -e "Failure reclibrating bam"
     exit 1
    fi

    mkdir 13_final_bam

    mv ./10_recalibrated_bam/recalibrated.bam ./13_final_bam/<%= sample_name %>.bam
    mv ./10_recalibrated_bam/recalibrated.bai ./13_final_bam/<%= sample_name %>.bam.bai
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

end


class AnalysisTemplaterApp
  VERSION       = "2.3.0"
  REVISION_DATE = "20130602"
  AUTHOR        = "Stuart Glenn <Stuart-Glenn@omrf.org>"
  COPYRIGHT     = "Copyright (c) 2012-2013 Oklahoma Medical Research Foundation"

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
    @default_config = (@config["DEFAULT"] || [{:run => nil, :lane => nil, :bwa_ref => nil, :gatk_ref => nil, :snp_rod => nil,:opts=>{:skip_btangs => false, :skip_dupes =>false, :skip_vcf => false, :skip_indel_realign => false, :reduce_reads=>false}}]).first

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
    end
    if @options.debug
      puts data.inspect
      #exit 0
      puts AnalysisTemplate.new(@default_config,sample_name,data)
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
      f.puts AnalysisTemplate.new(@default_config,sample_name,data)
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
    cmd = "qsub -o logs -sync y -b y -V -j y -cwd -q ngs.q -m e -N a_#{sample_name}_full ./analyze.sh"
    cmd = "./analyze.sh" if @options.run_local
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
  -l, --local            Run the analyze script locally, not with initial SGE submit
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

# It is easier to use variables for thse paths
GATK_REF=<%= @data.first[:gatk_ref] || @default_config[:gatk_ref] %>
GATK_DBSNP=<%= @data.first[:snp_rod] || @default_config[:snp_rod] || 'nil' %>
GATK_BIN=`which gatk`
GATK_BASE=`dirname ${GATK_BIN}`"/.."

SAMPLE="<%= @sample_name %>"

<%=
  fastq_file_list(@sample_name,@data)
%>

mkdir tmp
export GATK_JAVA_OPTS="-Djava.io.tmpdir=./tmp"

if [ ! -e 05_dup_marked/cleaned.bam ]; then
# initial PCR clean by 'bins' via btangs
<%=
  clean_commands(@sample_name,@data)
%>

# fastqc info
mkdir qc
qsub -l virtual_free=2G,h_vmem=4G -o logs -b y -V -j y -cwd -q ngs.q -N a_<%= @sample_name %>_qc fastqc -o qc <%= fastq_shell_vars() %>

# setup input sams, will get illumina scores to standard sanger
#mkdir 00_inputs
export JAVA_MEM_OPTS="-Xmx16G"

# Now take those & actually generate the alignment SAM output, paired or single
mkdir 03_sorted_bams
<%=
  bwa_alignment_command(@sample_name,@data)
%>

if [ "$?" -ne "0" ]; then
  echo -e "Failure with bwa alignment"
  exit 1
fi


# Going forward let us play with BAM files
# We will sort the individual aligned SAM files so we can merge, sort, & mark dupes in one action
#mkdir 03_sorted_bams
#qsub -l virtual_free=4G,h_vmem=56G -o logs -sync y -t 1-<%= total_number_input_sequenced_lanes() %> -b y -V -j y -cwd -q ngs.q -N a_<%= @sample_name %>_sort_sams sort_sam_qsub_tasked.rb 03_sorted_bams <%= input_sam_bam_files("02_bwa_alignment","sam") %>

#if [ "$?" -ne "0" ]; then
  #echo -e "Failure sorting aligned sams"
  #exit 1
#fi

rm -rf 00_inputs 01_bwa_aln_sai 02_bwa_alignment


# Now we might have had many input SAMs, so let us merge those all into a single BAM using picard
# While we do that we shall also sort it, mark possible duplicates & make an index
mkdir 05_dup_marked
<%= mark_dupes_or_skip(@sample_name,@data) %>

<%= (@data.first.has_key?(:keep_unaligned) && @data.first[:keep_unaligned]) ? "": "rm -rf 03_sorted_bams" %> 


# TODO stop here if pre_gatk only
if [ "$PRE_GATK_ONLY" == "Y" ]; then
rm -rf 00_inputs \
01_bwa_aln_sai \
02_bwa_alignment <%= (@data.first.has_key?(:keep_unaligned) && @data.first[:keep_unaligned]) ? "": "03_sorted_bams" %> \
tmp
<%=
  cleanup_cleaned_fastq_files(@sample_name)
%>

rm -f qc/*.zip

touch finished.txt
exit 0
fi
else
  rm -f finished.txt
  mkdir tmp
fi #if 05_dup_marked/cleaned.bam already existed
# start here if pre_gatk already done

<%= indel_realign(@sample_name,@data) %>

<%= covariate_or_final(@sample_name,@data) %>

# fastqc info
qsub -l virtual_free=2G,h_vmem=4G -o logs -b y -V -j y -cwd -q ngs.q -N a_<%= @sample_name %>_qc fastqc -o qc ./13_final_bam/<%= @sample_name %>.bam

<%= reduce_reads(@sample_name,@data) %>

<%= variant_call(@sample_name,@data) %>

# Clean up after ourselves
rm -rf 00_inputs \
01_bwa_aln_sai \
02_bwa_alignment <%= (@data.first.has_key?(:keep_unaligned) && @data.first[:keep_unaligned]) ? "": "03_sorted_bams" %> \
05_dup_marked \
06_intervals \
07_realigned_bam \
08_uncalibated_covariates \
10_recalibrated_bam \
11_calibated_covariates 


<%=
  if (@data.first.has_key?(:keep_unaligned) && @data.first[:keep_unaligned]) then
    cmd = "qsub -o logs -b y -V -j y -cwd -q ngs.q -m e -N a_#{@sample_name}_unaligned_extract ./extract_unaligned.sh"
    cmd
  end
%>

# gzip & clean up some of the cleaned input as we keep some it for now, but don't need it fullsized
if [ "$PRE_GATK_ONLY" != "Y" ]; then
<%=
  cleanup_cleaned_fastq_files(@sample_name)
%>
fi

rm -f qc/*.zip
rm -rf tmp

touch finished.txt
