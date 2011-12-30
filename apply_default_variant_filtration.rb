#!/usr/bin/env ruby1.9
#
# apply_default_variant_filtration.rb
# Created by Stuart Glenn on 2011-06-09T15:39:34-0500 
#
# == Synopsis
# Quick script to wrap the work of running our standard GATK variant filters
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
output_vcf = File.join(base_dir,"#{File.basename(file,extension)}-filtered#{extension}")

if File.exist?(output_vcf) then
  STDERR.puts "A file with the output name we would use (#{output_vcf}) already exists. I won't overwrite"
  exit 1
end

cmd = <<-EOF
gatk -et NO_ET -T VariantFiltration \\
-R #{reference} \\
-o #{output_vcf}  \\
-V #{input_vcf} \\
--filterExpression " MQ < 40.0 || QD < 2.0 || FS > 60.0 || HaplotypeScore > 13.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0"  \\
--filterName "GATKStandard"
EOF

puts cmd
exec(cmd)
