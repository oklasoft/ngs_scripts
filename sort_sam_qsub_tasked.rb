#!/usr/bin/env ruby1.8

output_base = ARGV.shift

index = (ENV['SGE_TASK_ID']||1).to_i - 1

data = ARGV[index]

output = File.join(output_base,"#{index}.bam")

cmd = "picard SortSam TMP_DIR=./tmp VALIDATION_STRINGENCY=LENIENT MAX_RECORDS_IN_RAM=3000000 OUTPUT=#{output} SORT_ORDER=coordinate INPUT=#{data}"

puts cmd
STDOUT.flush

exec cmd
