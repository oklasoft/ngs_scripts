#!/usr/bin/env ruby1.8

output_base = ARGV.shift
reference = ARGV.shift
index = ENV['SGE_TASK_ID'].to_i - 1
slots = ENV['NSLOTS']

inputs = ARGV.each_slice(2).to_a[index]

output = File.join(output_base,"#{File.basename(inputs.last,".bam")}-#{inputs.first}.sai")

cmd = "bwa aln -t #{slots} -f #{output} #{reference} -b#{inputs.first} #{inputs.last}"

puts cmd
STDOUT.flush

exec cmd
