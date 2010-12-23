#!/usr/bin/env ruby1.8

args = ARGV.clone
# args.shift
output_base = args.shift
platform = args.shift
sample = args.shift

index = (ENV['SGE_TASK_ID']||1).to_i - 1


tag = args[index]

cmd = "add_read_group_to_sam.pl -i #{index}.sam -R -r #{tag} -s #{sample} -p #{platform}"

puts cmd
STDOUT.flush

#exec cmd
