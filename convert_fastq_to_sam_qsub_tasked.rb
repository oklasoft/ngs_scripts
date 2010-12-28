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
  data[:id] = args.shift
  data[:unit] = args.shift
  
  data[:inputs] = ["#{args.shift}"]  
  data[:inputs] << ["#{args.shift}"] if "paired" == data[:mode]
  
  groups << data
end

data = groups[index]

output = File.join(output_base,"#{index}.sam")

cmd = "picard FastqToSam VALIDATION_STRINGENCY=LENIENT MAX_RECORDS_IN_RAM=1000000 QUALITY_FORMAT=Illumina OUTPUT=#{output} READ_GROUP_NAME=#{data[:id]} SAMPLE_NAME=#{sample} PLATFORM=#{platform} PLATFORM_UNIT=#{data[:unit]} FASTQ=#{data[:inputs].shift}"

unless data[:inputs].empty? then
  cmd += " FASTQ2=#{data[:inputs].shift}"
end

puts cmd
STDOUT.flush

unless system(cmd)
  STDERR.puts "Failed with picard"
  exit 1
end

output_bam = File.join(output_base,"#{index}.bam")

cmd = "samtools view -bhS #{output} > #{output_bam}"

puts cmd
STDOUT.flush

exec cmd
