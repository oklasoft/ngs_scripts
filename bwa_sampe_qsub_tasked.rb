#!/usr/bin/env ruby1.8

output_base = ARGV.shift
reference = ARGV.shift
index = (ENV['SGE_TASK_ID']||1).to_i - 1

inputs = ARGV.each_slice(5).to_a[index]
tag = inputs.shift
output = File.join(output_base,"#{File.basename(inputs.first,"_1.sai")}.sam")

cmd = "bwa sampe -r \"#{tag}\" -f #{output} #{reference} #{inputs.join(" ")}"

puts cmd
STDOUT.flush

exec cmd
