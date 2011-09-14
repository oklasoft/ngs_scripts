#!/usr/bin/env ruby
#
# comparae_vcf_to_plink.rb
# Created by Stuart Glenn on 20110914
#
# == Synopsis
# Quick script wrapping up the steps for taking VCF & --genome compare to some plink data
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

require 'ostruct'
require 'optparse'
require 'erb'
require 'yaml'
require 'open4'
require 'fileutils'

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
  # exit(status.exitstatus)
end

def output_version(out)
  out.puts "#{File.basename($0)} 1.0"
end

def output_help(out)
  out.puts "#{File.basename($0)} -p PLINK_BASE_FILE -s SNPS.TXT -o OUTPUT_BASE_DIR -i INPUT.VCF"
  out.puts ""
  out.puts <<-EOF
  -v, --version             Print version info
  -h, --help                Print this help
  -V, --verbose             Enable verbosity
  -D, --debug               Enable debuging
  -o, --output DIR          Save output in DIR
  -i, --input VCF           Input VCF to be checked
  -p, --plink BASE_PED      Base path/name to plink ped/map for merging
  -s, --snps FILE.TXT       List of SNP names that will be used for check
EOF
end

def parse_opts(args)
  @options = OpenStruct.new(
    :output_base_dir => Dir.pwd,
    :verbose => false,
    :debug => false,
    :input_vcf_path => nil,
    :plink_data_path => nil,
    :snp_list_path => nil
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

    o.on("-p","--plink", "=REQUIRED") do |path|
      @options.plink_data_path = File.expand_path(path)
    end

    o.on("-s","--snp", "=REQUIRED") do |path|
      @options.snp_list_path = File.expand_path(path)
    end
  end

  opts.parse!(args)  

  unless @options.input_vcf_path && @options.plink_data_path && @options.snp_list_path
    $stderr.puts "Need at least an input VCF, snp list file & plink data"
    $stderr.puts ""
    output_help($stderr)
    exit 1
  end  
end

def has_plink?
  result = ""
  status = Open4::open4("which","plink") do |pid, stdin, stdout, stderr|
    stdin.close
    result = stdout.read.strip unless stdout.eof?
  end
  if 0 == status.exitstatus
    return result
  else
    return false
  end
end

parse_opts(ARGV)

unless has_plink?
  $stderr.puts "Can't find plink"
  exit(1)
end

base_vcf_name = File.basename(File.basename(@options.input_vcf_path,".gz"),".vcf")

Dir.chdir(@options.output_base_dir) do
  puts "Comparison Starting #{$$}" if @options.verbose
  work_dir = "#{Time.now.strftime("%Y%m%d")}_compare_vcf_to_plink_work_dir.#{$$}"
  raise "Can't make work dir" unless Dir.mkdir(work_dir)
  Dir.chdir(work_dir) do
    get_copy_of_vcf
    uncompress_vcf
    convert_vcf_to_plink
    extract_subset_of_snp_from_plink
    merge_compare_set_to_vcf_set
    flip_strand_vcf_set
    merge_compare_set_to_vcf_set
    exclude_any_remaining_problem_snps
    merge_compare_set_to_vcf_set
    genome_compare
    report_results
  end
end
