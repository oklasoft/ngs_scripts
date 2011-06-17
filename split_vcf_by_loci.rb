#!/usr/bin/env ruby1.9 -w -Ku
#
# split_vcf_by_loci.rb
# Created by Stuart Glenn on 2011-06-13T13:39:57-0500 
#
# == Synopsis
# Quick script to wrap the work of splitting a VCF into sub VCFs by loci
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

def run_command(cmd,msg="Command failed")
  puts cmd
  unless system(cmd)
    STDERR.puts "#{msg}: #{$?}"
    exit $?.exitstatus
  end
end

def subset_vcf_to_vcf_around(input_vcf,output_vcf,chr,start,stop)
  cmd = <<-EOF
  vcftools --vcf #{input_vcf} \\
  --out #{output_vcf} \\
  --chr #{chr} \\
  --from-bp #{start} \\
  --to-bp #{stop} \\
  --remove-filtered GATKStandard \\
  --remove-filtered HARD_TO_VALIDATE \\
  --remove-filtered LowQual \\
  --remove-filtered SnpCluster \\
  --keep-INFO-all \\
  --recode
EOF
  run_command(cmd)
end

def get_input_file(arg,name)
  file = arg
  unless file then
    STDERR.puts "Missing input #{name} file"
    STDERR.puts "Usage: #{__FILE__} <INPUT_VCF> <FILE_OF_LOCI> [BASE_OUTPUT_DIR]"
    exit 1
  end

  file = File.expand_path(file)
  unless File.readable?(file)
    STDERR.puts "The input #{name} file is not actually readable"
    exit 1
  end  
  return file
end


def main()
  input_vcf = get_input_file(ARGV.shift,"VCF")
  loci_file = get_input_file(ARGV.shift,"Loci")
  output_base = get_input_file(ARGV.shift || Dir.pwd,"Output base")

  (base_dir,file) = File.split(input_vcf)
  extension = File.extname(file)
  if ".gz" == extension
    extension = "#{File.extname(File.basename(file,extension))}"
  end

  IO.foreach(loci_file) do |locus_line|
    next if locus_line =~ /^#/ || locus_line =~ /^$/
    (chr,start,stop,name) = locus_line.chomp.split(/\s+/)
    STDERR.puts "Splitting some stuff out for #{name} in #{chr}:#{start}-#{stop}"
    
    output_vcf = File.join(output_base,name,"#{File.basename(file,extension)}-#{name}")
    if File.exist?(output_vcf) then
      STDERR.puts "A file with the output name we would use (#{output_vcf}) already exists. I won't overwrite"
      exit 1
    end

    unless Dir.mkdir(File.join(output_base,name))
      STDERR.puts "Failed to make the output dir for #{name} in #{output_base}"
    end
    
    subset_vcf_to_vcf_around(input_vcf,output_vcf,chr,start,stop)
  end
end

main()