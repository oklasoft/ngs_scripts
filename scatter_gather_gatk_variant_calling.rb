#!/usr/bin/env ruby -w -Ku
#
# analyze_sequence_to_snps.rb
# Created by Stuart Glenn on 2011-06-07
#
# == Synopsis
# Quick script to wrap the work of running GATK scatter/gathter type || UnifiedGenotyper
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

@options = OpenStruct.new(
  :output_base => Dir.pwd,
  :scatter_limit => 100,
  :verbose => true,
  :bam_list => nil,
  :interval_list => nil
)

@options.bam_list = ARGV.shift
@options.interval_list = ARGV.shift
@options.output_base = ARGV.shift unless ARGV.empty?
@options.scatter_limit = ARGV.shift.to_i unless ARGV.empty?

unless @options.bam_list && @options.interval_list
  STDERR.puts "Need to at least give us a file listing bams & the master intervals"
  exit 1
end



STDERR.puts "#{@options}" if @options.verbose
name_base = File.split(@options.bam_list)[-1]

Dir.chdir(@options.output_base) do
  puts "Scatter/Gather Variants Starting #{$$}"
  work_dir = "#{Time.now.strftime("%Y%m%d")}_scatter_gather_work.#{$$}"
  raise "Can't make work dir" unless Dir.mkdir(work_dir)
  to_joins = []
  Dir.chdir(work_dir) do
    raise "Can't make log dir" unless Dir.mkdir("logs")
    lines = IO.readlines(@options.interval_list).map {|l| l.chomp}
    intervals_per_scatter = (lines.size/@options.scatter_limit.to_f).ceil.to_i
    slice = 0
    lines.each_slice(intervals_per_scatter) do |sliced_intervals|
      slice = slice.to_s.rjust(@options.scatter_limit.to_s.length,"0")
      
      STDERR.puts "#{slice}\t#{sliced_intervals.join("; ")}" if @options.verbose

      sliced_interval_file = "sliced_interval_#{slice}.interval_list"
      File.open(sliced_interval_file,"w") do |f|
        f.puts sliced_intervals.join("\n")
      end

      output_file = File.join(Dir.pwd,"#{name_base}_#{slice}_variants.vcf")
      input_file = File.join(Dir.pwd,sliced_interval_file)

      cmd = "qsub -m e -pe threaded 6 -o logs -b y -V -j y -cwd -q all.q -N #{name_base}_variants_#{slice} \
gatk -et NO_ET -T UnifiedGenotyper -glm BOTH -nt 6 \
-R /Volumes/hts_core/Shared/homo_sapiens_36.1/chr_fixed/hg18.fasta \
-I #{@options.bam_list} \
-o #{output_file} \
-D /Volumes/hts_core/Shared/dbsnp/dbsnp_130_hg18.rod \
-L #{input_file}"
      puts cmd
      system cmd
      sleep(10)
      slice = slice.to_i + 1
      
      to_joins << "-B:#{input_file},VCF #{output_file}"
    end #each slice
  end #work_dir

  cmd = "qsub -m e -o logs -b y -V -j y -cwd -q all.q -N #{name_base}_variants_merge \
gatk -et NO_ET -T CombineVariants \
-variantMergeOptions UNION \
-genotypeMergeOptions UNSORTED \
-R /Volumes/hts_core/Shared/homo_sapiens_36.1/chr_fixed/hg18.fasta \
-o #{name_base}_variants.vcf \
#{to_joins.join(" ")}
"
  File.open("#{work_dir}-joiner.sh","w") do |f|
    f.puts "#!/usr/bin/env bash"
    f.puts
    f.puts "module load sge"
    f.puts "module load gatk"
    f.puts cmd
  end
  
end #base_dir
