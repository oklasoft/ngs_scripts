#!/usr/bin/env ruby1.8

output_base = ARGV.shift
reference = ARGV.shift
index = (ENV['SGE_TASK_ID']||1).to_i - 1

input = ARGV[index]

output = File.join(output_base,"#{File.basename(input,".sam")}.bam")
sorted_output = File.join(output_base,"#{File.basename(input,".sam")}-sorted")

cmd = "samtools view -bt #{reference} -i #{output} #{input}"

puts cmd
unless system(cmd)
  STDERR.puts "Failed with samtools import"
  exit 1
end
cmd = "samtools sort #{output} #{sorted_output}"

puts cmd
unless system(cmd)
  STDERR.puts "Failed with samtools sort"
  exit 1
end
cmd = "rm #{output} && mv #{sorted_output}.bam #{output}"
puts cmd
STDOUT.flush

exec cmd
