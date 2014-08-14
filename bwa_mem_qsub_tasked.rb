#!/usr/bin/env ruby1.9

require 'tmpdir'

args = ARGV.clone
# args.shift
tmp_base = args.shift
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

return_val = -1
Dir.mktmpdir(index.to_s,tmp_base) do |tmp_prefix_dir|
  cmd = "bwa mem -v 1 -M -t #{threads.to_i-2} -R \"#{tag}\" #{reference} #{data[:inputs].join(" ")}"
  cmd += "|samtools view -Shu - | samtools sort -@ 2 -m 4G -o - #{tmp_prefix_dir}/#{index} > #{output}"

  puts cmd
  STDOUT.flush
  if system "/bin/bash", "-o", "pipefail", "-o", "errexit", "-c", cmd
    return_val = 0
  else
    return_val = $?.exitstatus
  end
end

exit return_val
