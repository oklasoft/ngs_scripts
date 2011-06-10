#!/usr/bin/env ruby1.9 -w -Ku
#
# phase_vcf_with_beagle.rb
# Created by Stuart Glenn on 2011-06-10T09:18:14-0500  
#
# == Synopsis
# Quick script to wrap the work of running our beagle phasing
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

BEAGLE_LIKELIHOOD_INPUT = "beagle_likelihood_input"
BEAGLE_OUTPUT_BASENAME = "beagle_phased_output"

def run_command(cmd,msg="Command failed")
  puts cmd
  unless system(cmd)
    STDERR.puts "#{msg}: #{$?}"
    exit $?.exitstatus
  end
end

def produce_likelihood_input_for_beagle(input_vcf,reference)
  cmd = <<-EOF
  gatk -et NO_ET -T ProduceBeagleInput \\
  -R #{reference} \\
  -o #{BEAGLE_LIKELIHOOD_INPUT}  \\
  -B:variant,VCF #{input_vcf}
EOF
  run_command(cmd,"Producing Likelihood for Beagle Failed")
end

def phase_with_beagle()
  cmd = <<-EOF
  beagle.sh \\
  like=#{BEAGLE_LIKELIHOOD_INPUT}  \\
  -out=#{BEAGLE_OUTPUT_BASENAME}
EOF
  run_command(cmd,"Phasing with Beagle Failed")
end

def uncompress_beagle_ouput()
  %w{phased gprobs}.each do |f|
    cmd="gunzip #{BEAGLE_OUTPUT_BASENAME}.beagle_output.#{f}.gz"
    run_command(cmd,"Failed uncompressing #{f}")
  end
end


def convert_beagle_to_vcf(input_vcf,reference,output_vcf)
  cmd = <<-EOF
  gatk -T BeagleOutputToVCF \\
  -R #{reference} \\
  -o #{output_vcf} \\
  -B:variant,VCF #{input_vcf} \\
  -B:beagleR2,BEAGLE #{BEAGLE_OUTPUT_BASENAME}.beagle_output.r2 \\
  -B:beaglePhased,BEAGLE #{BEAGLE_OUTPUT_BASENAME}.beagle_output.phased \\
  -B:beagleProbs,BEAGLE #{BEAGLE_OUTPUT_BASENAME}.beagle_output.gprobs
EOF
  run_command(cmd,"Failed converting beagle back to VCF")
end

def main()
  input_vcf = ARGV.shift
  reference = ARGV.shift || "/Volumes/hts_core/Shared/homo_sapiens_36.1/chr_fixed/hg18.fasta"

  unless input_vcf then
    STDERR.puts "Missing input VCF file"
    exit 1
  end

  input_vcf = File.expand_path(input_vcf)

  unless File.readable?(input_vcf)
    STDERR.puts "The input VCF file is not actually readable"
    exit 1
  end

  (base_dir,file) = File.split(input_vcf)
  extension = File.extname(file)
  if ".gz" == extension
    extension = "#{File.extname(File.basename(file,extension))}"
  end
  output_vcf = File.join(base_dir,"#{File.basename(file,extension)}-phased#{extension}")

  if File.exist?(output_vcf) then
    STDERR.puts "A file with the output name we would use (#{output_vcf}) already exists. I won't overwrite"
    exit 1
  end

  work_dir = File.join(base_dir,"#{Time.now.strftime("%Y%m%d")}_scatter_gather_work.#{$$}")
  raise "Can't make work dir" unless Dir.mkdir(work_dir)
  Dir.chdir(work_dir) do
    produce_likelihood_input_for_beagle(input_vcf,reference)
    phase_with_beagle()
    uncompress_beagle_ouput()
    convert_beagle_to_vcf(input_vcf,reference,output_vcf)
  end
end

main()