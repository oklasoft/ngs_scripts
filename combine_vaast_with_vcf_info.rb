#!/usr/bin/env ruby1.9
#
# combine_vaast_with_vcf_info.rb
# Created by Stuart Glenn on 2011-11-15
#
# == Synopsis
# Script to combine the vaast output results with some VCF info data
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

$:.unshift(File.join(File.dirname(__FILE__),"..","lib"))
require "vaast_database"

vcf_info_file = ARGV.shift
vaast_files = ARGV.map {|f| File.expand_path(f)}

unless vcf_info_file && vaast_files.size >= 1
  $stderr.puts "Missing input files:"
  $stderr.puts "Usage:"
  $stderr.puts "#{File.basename(__FILE__)} VCF_INFO_FILE VAAST_FILE(S)"
  exit(1)
end

vaast = OMRF::VaastDatabaseOuputFiles.new(vaast_files)

type_priorities = [:tu, :tr, :br]

IO.foreach(vcf_info_file) do |line|
  if 1 == $.
    print "#{line.chomp}\tvaast_nm_id\tvaast_gene_name\tvaast_score\tvaast_rank\tvaast_genome_permutation_p"
    print "\tvaast_genome_permutation_0.95_lower\tvaast_genome_permutation_0.95_upper"
    puts "\tvaast_type\tvaast_score\tvaast_change"
  else
    parts = line.chomp.split(/\t/)
    chrom = parts[0]
    position = parts[1]
    features = vaast.find(chrom,position.to_i)
    parts += [features[:id],features[:gene_name],features[:score],features[:rank],
      features[:genome_permutation_p],features[:genome_permutation_95_lower],
      features[:genome_permutation_95_upper]]

    types = features[:types] || {}
    type_name = nil
    type_priorities.each do |t|
      if types[t]
        type_name = t
        break
      end
    end
    type = types[type_name] || {}
    parts += [type_name, type[:score], type[:change]]
    puts parts.join("\t")
  end
end
