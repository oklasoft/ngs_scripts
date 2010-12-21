#!/usr/bin/env ruby1.8

output_base = ARGV.shift
reference = ARGV.shift
index = ENV['SGE_TASK_ID'].to_i - 1

input = ARGV[index]

output = File.join(output_base,"#{File.basename(input,".fastq")}.sai")

cmd = "bwa aln -t 12 -f #{output} #{reference} #{input}"

puts cmd
STDOUT.flush

exec cmd
