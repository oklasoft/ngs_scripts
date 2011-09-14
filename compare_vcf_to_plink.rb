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

def fail(msg)
  $stderr.puts msg
  exit(1)
end

def get_copy_of_vcf()
  FileUtils.cp(@options.input_vcf_path,File.basename(@options.input_vcf_path))
  File.readable?(File.basename(@options.input_vcf_path)) #why you no have return value cp
end

def uncompress_vcf()
  cmd = ["gunzip",File.basename(@options.input_vcf_path)]
  run_command("uncomress VCF",cmd)
end

def convert_vcf_to_plink()
  cmd = ["vcftools","--plink","--vcf","#{@options.base_name}.vcf","--out",@options.base_name]
  run_command("vcf to plink",cmd)
end

def extracted_name
  "#{@options.base_name}_extracted"
end

def flipped_name
  "#{@options.base_name}_extracted_flipped"
end

def excluded_name
  "#{@options.base_name}_extracted_flipped_excluded"
end

def extract_subset_of_snp_from_plink()
  cmd = ["plink","--file",@options.base_name,"--extract",@options.snp_list_path,"--recode","--out",extracted_name]
  run_command("extract SNPs",cmd)
end

def merge_data_set_with_as(a,b,output)
  cmd = [
    "plink",
    "--file",
    a,
    "--merge",
    "#{b}.ped",
    "#{b}.map",
    "--recode",
    "--out",
    output
  ]
  run_command("plink merge",cmd)
end

def flip_strand_for_in(snp_list,plink_set)
  cmd = ["plink","--file",plink_set,"--flip",snp_list,"--recode","--out",flipped_name()]
  run_command("plink strand flip",cmd)
end

def exclude_from_in(snp_list,plink_set)
  File.open("bad_snps_to_extract.txt","w") do |f|
    IO.foreach(snp_list) do |line|
      parts = line.chomp.split(/\t/)
      f.puts parts[1]
    end
  end
  cmd = ["plink","--file",plink_set,"--extract","bad_snps_to_extract.txt","--recode","--out",excluded_name()]
  run_command("plink strand flip",cmd)
end

def merge_data_sets_fixing_flips()
  unless merge_data_set_with_as(@options.plink_data_path,extracted_name(),"merge")
    # a failed first merge means some strand flipping is needed most likely
    if File.size?("merge.missnp")
      flip_strand_for_in("merge.missnp",extracted_name()) || fail("Unable to flip strands")
      File.unlink("merge.missnp")
    else
      fail("First merge failed & there were no strands to flip")
    end
    unless merge_compare_set_to_vcf_set(@options.plink_data_path,flipped_name(),"merge")
      # second merge failure is likely due then to extra screw snps, exclude those suckers
      if File.size?("merge.missnp")
        exclude_from_in("merge.missnp",flipped_name()) || fail("Unable to exclude problem SNPs")
        merge_data_set_with_as(@optopns.plink_data_path,excluded_name(),"merge")
      else
        fail("Second merge failed & there were no SNPs to extract")
      end
    end
  end
end

def run_comparison()
  get_copy_of_vcf() || fail("Could not get VCF copy")
  uncompress_vcf || fail("Could not uncompress VCF")
  convert_vcf_to_plink || fail("Could not convert VCF to plink")
  extract_subset_of_snp_from_plink || fail("Could not extract SNP subset")
  merge_data_sets_fixing_flips() || fail("Could not merge")
  genome_compare
  report_results
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

def input_files_readable?
  {
    @options.input_vcf_path  => "input VCF",
    @options.snp_list_path  => "SNP list",
    "#{@options.plink_data_path}.ped" => "plink ped file",
    "#{@options.plink_data_path}.map" => "plink map file"
  }.each do |file,msg|
    unless File.readable?(file)
      $stderr.puts "Sorry I can't read the #{msg}"
      return false
    end
  end
  return true
end

def has_cmd?(cmd)
  result = ""
  status = Open4::open4("which",cmd) do |pid, stdin, stdout, stderr|
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

unless input_files_readable? 
  $stderr.puts "Can't work without input"
  exit(1)
end

%w/plink vcftools gunzip/.each do |cmd|
  unless has_cmd?(cmd)
    $stderr.puts "Can't find #{cmd}"
    exit(1)
  end
end

@options.base_name = File.basename(File.basename(@options.input_vcf_path,".gz"),".vcf")

Dir.chdir(@options.output_base_dir) do
  puts "Comparison Starting #{$$}" if @options.verbose
  work_dir = "#{Time.now.strftime("%Y%m%d")}_compare_vcf_to_plink_work_dir.#{$$}"
  raise "Can't make work dir" unless Dir.mkdir(work_dir)
  Dir.chdir(work_dir) do
    run_comparison()
  end
end
