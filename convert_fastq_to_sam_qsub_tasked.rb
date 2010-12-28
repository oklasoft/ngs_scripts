#!/usr/bin/env ruby1.8

args = ARGV.clone
# args.shift
output_base = args.shift
platform = args.shift
sample = args.shift

index = (ENV['SGE_TASK_ID']||1).to_i - 1

groups = []

while a = args.shift
  data = {}
  if "paired" == a then
    data[:mode] = "paired"
  else
    data[:mode] = "single"
  end
  data[:unit] = args.shift
  data[:id] = args.shift
  
  data[:inputs] = ["#{args.shift}"]  
  data[:inputs] << ["#{args.shift}","#{args.shift}"] if "paired" == data[:mode]
  
  groups << data
end

data = groups[index]

output = File.join(output_base,"#{index}.sam")

cmd = "picard FastqToSam QUALITY_FORMAT=Illumina OUTPUT=#{output} READ_GROUP_NAME=#{data[:id]} SAMPLE_NAME=#{sample} PLATFORM=#{platform} PLATFORM_UNIT=#{dat[:unit]} FASTQ=#{data[:inputs].first}"

unless data[:inputs].empty? then
  cmd += " FASTQ2=#{data[:inputs].shift}"
end

puts cmd
STDOUT.flush

exec cmd
