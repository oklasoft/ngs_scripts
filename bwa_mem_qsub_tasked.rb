#!/usr/bin/env ruby1.9

args = ARGV.clone
# args.shift
output_base = args.shift
reference = args.shift
index = (ENV['SGE_TASK_ID']||1).to_i - 1
threads = (ENV['NSLOTS']) || 1

groups = [] #ARGV.each_slice(5).to_a[index]

while a = args.shift
  data = {}
  if "paired" == a then
    data[:mode] = "paired"
  else
    data[:mode] = "single"
  end
  data[:tag] = args.shift
  data[:inputs] = ["#{args.shift}"]
  
  data[:inputs] << ["#{args.shift}"] if "paired" == data[:mode]
  
  groups << data
end

data = groups[index]

tag = data[:tag]

output = File.join(output_base,"#{index}.bam")

cmd = "bwa mem -M -t #{threads} -R \"#{tag}\" #{reference} #{data[:inputs].join(" ")}"
cmd += "| picard SortSam TMP_DIR=./tmp VALIDATION_STRINGENCY=LENIENT MAX_RECORDS_IN_RAM=3000000 OUTPUT=#{output} SORT_ORDER=coordinate INPUT=/dev/stdin"

puts cmd
STDOUT.flush
#STDOUT.reopen(File.open(output, 'w'))
#exec cmd
if system cmd
  exit 0
else
  exit $?
end
