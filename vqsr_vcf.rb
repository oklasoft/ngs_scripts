#!/usr/bin/env ruby
#
# vqsr_vcf.rb
# Created by Stuart Glenn on 2011-08-05
#
# == Synopsis
# Quick script wrapping up the steps for VQSR in our pipeline
#
# ==Author
# Stuart Glenn <Stuart-Glenn@omrf.org>
#
# ==Copyright
#  Copyright (c) 2011,2016 Stuart Glenn, Oklahoma Medical Research Foundation. (OMRF)
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

require 'ostruct'
require 'optparse'
require 'erb'
require 'yaml'
require 'open4'
require 'fileutils'

# format out the -B options for gatk based on the vqsr data
def vqsr_training_data_gatk_opts(vqsr_data)
  opts = []
  vqsr_data[:training_sets].each do |set|
    opts << "-resource:#{set[:name]},#{set[:type]},#{set[:params]}"
    opts << set[:path]
  end
  vqsr_data[:annotations].each do |an|
    opts << "-an"
    opts << an
  end
  opts
end

def apply_recalibration_extra_opts(vqsr_data)
  extra_opts_for_key(vqsr_data,:apply_recalibration_opts)
end

def variant_recalibrator_extra_opts(vqsr_data)
  extra_opts_for_key(vqsr_data,:variant_recalibrator_opts)
end

def extra_opts_for_key(vqsr_data,key)
  opts = []
  return opts unless vqsr_data[key]
  vqsr_data[key].each do |opt,arg|
    opts << "--#{opt}"
    opts << arg.to_s
  end
  opts
end

def vqsr_output_file(base_vcf_name,type)
  File.join(File.expand_path("."),"#{base_vcf_name}.#{type}")
end

def run_command(name,cmd)
  opts = {
    'stdin' => nil, 
    'stdout' => $stdout, 
    'stderr' => $stderr, 
    'raise' => false, 
    'quiet' => true
  }
  puts "Running: #{cmd.join(" ")}" if @options.verbose
  pre = Hash.new
  pre[:stdout] = $stdout.sync
  pre[:stderr] = $stderr.sync
  $stdout.sync = true
  $stderr.sync = true
  status = Open4::spawn(cmd, opts)
  if nil == status
    $stderr.puts "Unable to start #{name}"
    $stderr.flush
    return false
  end
  $stdout.sync = pre[:stdout]
  $stderr.sync = pre[:stderr]
  return true if 0 == status.exitstatus
  $stderr.puts "#{name} exited with error: #{status.inspect}"
  $stderr.flush
  return false
  exit(status.exitstatus)
end

def variant_recalibrator(data,input_vcf,base_vcf_name)
  cmd = %w/gatk -T VariantRecalibrator -R/ + [@config[:gatk_ref]]
  cmd += vqsr_training_data_gatk_opts(data)
  cmd += variant_recalibrator_extra_opts(data)
  cmd += ["-input", input_vcf]
  cmd += ["-recalFile", vqsr_output_file(base_vcf_name,"recal")]
  cmd += ["-tranchesFile", vqsr_output_file(base_vcf_name,"tranches")]
  cmd += ["-rscriptFile", vqsr_output_file(base_vcf_name,"R")]
  run_command("VariantRecalibrator",cmd)
end

def make_plots(data,input_vcf,base_vcf_name)
  cmd = %w/R CMD BATCH --no-save/
  cmd += [vqsr_output_file(base_vcf_name,"R"), vqsr_output_file(base_vcf_name,"Rout")]
  run_command("R",cmd)
end


def apply_recalibration(data,input_vcf,base_vcf_name)
  cmd = %w/gatk -T ApplyRecalibration -R/ + [@config[:gatk_ref]]
  cmd += apply_recalibration_extra_opts(data)
  cmd += ["-input", input_vcf]
  cmd += ["-recalFile", vqsr_output_file(base_vcf_name,"recal")]
  cmd += ["-tranchesFile", vqsr_output_file(base_vcf_name,"tranches")]
  cmd += ["-o", File.join(@options.output_base_dir,"#{base_vcf_name}-vqsr_recalibrated.vcf")]
  run_command("ApplyRecalibration",cmd)
end

# check if we have the options for doing vsqr on the VCF,
# we'll have some yaml ala:
# ---
# :vqsr_$MODE:
#   :training_sets:
#   - :name: hapmap
#     :type: VCF
#     :path: hapmap_3.3.b37.sites.vcf
#     :params: known=false,training=true,truth=true,prior=15.0
#   - :name: omni
#     :type: VCF
#     :path: 1000G_omni2.5.b37.sites.vcf
#     :params: known=false,training=true,truth=false,prior=12.0
#   :annotations:
#   - HaplotypeScore
#   - MQRankSum
#   - ReadPosRankSum
#   - FS
#   - MQ
def do_vqsr?(vqsr_data)
  if vqsr_data
    if vqsr_data[:training_sets] && vqsr_data[:training_sets].size > 0 &&
      vqsr_data[:annotations] && vqsr_data[:annotations].size >0
      return true
    end
  end
  return false
end

def load_config_data(config_path)
  begin
    config = YAML::load(File.open(config_path))
    @config = config['DEFAULT'].first || config.first.first
    $stderr.puts @config if @options.verbose
  rescue => err
    $stderr.puts "Error loading config file: #{err}"
    exit(1)
  end
end

def output_version(out)
  out.puts "#{File.basename(__FILE__)} Version: 4.5.2-clia Released: 20171029"
end

def output_help(out)
  out.puts "#{File.basename($0)} -c CONFIG_YAML -o OUTPUT_BASE_DIR -i INPUT.VCF"
  out.puts ""
  out.puts <<-EOF
  -v, --version             Print version info
  -h, --help                Print this help
  -V, --verbose             Enable verbosity
  -D, --debug               Enable debuging
  -o, --output DIR          Save output in DIR
  -i, --input VCF           Input VCF to be VQSRed
  -m, --mode SNP|INDEL      Run as either SNP or INDEL (or some other mode, used to find vqsr_opts)
  -c, --config YAML         A configuration YAML, has more specific flags for VQSR
EOF
end

def parse_opts(args)
  @options = OpenStruct.new(
    :output_base_dir => Dir.pwd,
    :verbose => false,
    :debug => false,
    :input_vcf_path => nil,
    :mode => nil,
    :config_path => nil

  )
  
  opts = OptionParser.new() do |o|
    o.on('-v','--version') { output_version($stdout); exit(0) }
    o.on('-h','--help') { output_help($stdout); exit(0) }
    o.on('-V', '--verbose')    { @options.verbose = true }
    o.on('-D', '--debug')    { @options.debug = true }

    o.on("-o","--output", "=REQUIRED") do |output|
      @options.output_base_dir = File.expand_path(output)
    end

    o.on("-i","--input", "=REQUIRED") do |input|
      @options.input_vcf_path = File.expand_path(input)
    end

    o.on("-c","--config", "=REQUIRED") do |path|
      @options.config_path = File.expand_path(path)
    end

    o.on("-m","--mode", "=REQUIRED") do |mode|
      @options.mode = mode
    end
  end

  opts.parse!(args)  

  unless @options.input_vcf_path && @options.config_path && @options.mode
    $stderr.puts "Need at least an input VCF & an input config file & mode"
    $stderr.puts ""
    output_help($stderr)
    exit 1
  end  
end

def gatk_path
  result = ""
  status = Open4::open4("which","gatk") do |pid, stdin, stdout, stderr|
    stdin.close
    result = stdout.read.strip unless stdout.eof?
  end
  if 0 == status.exitstatus
    return result
  else
    $stderr.puts "Can't find gatk"
    exit(1)
  end
end

parse_opts(ARGV)
load_config_data(@options.config_path)
@vqsr_opts = @config["vqsr_#{@options.mode}".to_sym]
unless do_vqsr?(@vqsr_opts)
  $stderr.puts "Missing VQSR config options in YAML file for #{@options.mode}"
  exit(1)
end
@gatk_base = File.dirname(File.dirname(gatk_path))
base_vcf_name = File.basename(File.basename(@options.input_vcf_path,".gz"),".vcf")

Dir.chdir(@options.output_base_dir) do
  puts "VQSR Starting #{$$}" if @options.verbose
  work_dir = "#{Time.now.strftime("%Y%m%d")}_vqsr_vcf_work.#{$$}"
  raise "Can't make work dir" unless Dir.mkdir(work_dir)
  Dir.chdir(work_dir) do
    unless variant_recalibrator(@vqsr_opts,@options.input_vcf_path,base_vcf_name) &&
      make_plots(@vqsr_opts,@options.input_vcf_path,base_vcf_name) &&
      apply_recalibration(@vqsr_opts,@options.input_vcf_path,base_vcf_name)
      exit(1)
    else
      FileUtils.move(vqsr_output_file(base_vcf_name,"R.pdf"),File.join("..","#{base_vcf_name}-vqsr_recalibrated.pdf"))
      if File.exists?(File.join(vqsr_output_file(base_vcf_name,"tranches.pdf")))
        FileUtils.move(vqsr_output_file(base_vcf_name,"tranches.pdf"),File.join("..","#{base_vcf_name}-vqsr_recalibrated-tranches.pdf"))
      end
    end
  end
end
