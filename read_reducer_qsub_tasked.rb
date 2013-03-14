#!/usr/bin/env ruby1.8

HUMAN_CHRS = %w/1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y MT/

ref = ARGV.shift
output_base = ARGV.shift
output_prefix = ARGV.shift
input_bam = ARGV.shift

unless output_base && output_prefix && input_bam && ref
  raise "Missing some option <ref> <output_base_dir> <output_bam_prefix> <input_bam_file>"
end

Dir.chdir(output_base)

cmd = %w/gatk -T ReduceReads -R/ + [ref,"-I",input_bam]

if ENV['SGE_TASK_ID']
  # we will do it by chr
  begin
    Dir.mkdir('by_chr')
  rescue
  end
  index = (ENV['SGE_TASK_ID']).to_i - 1
  chr = HUMAN_CHRS[index]
  output_prefix = File.join("by_chr","#{output_prefix}-#{chr}")
  cmd += ["-L",chr]
else
  cmd += ["-compress","9"]
end
cmd += ["-o","#{output_prefix}.bam"]

puts cmd.join(" ")
STDOUT.flush

exec *cmd
