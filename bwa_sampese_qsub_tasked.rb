#!/usr/bin/env ruby1.8

output_base = ARGV.shift
reference = ARGV.shift
index = (ENV['SGE_TASK_ID']||1).to_i - 1

groups = [] #ARGV.each_slice(5).to_a[index]

ARGV.each do |a|
  data << {}

  if "paired" == a then
    data[:mode] = "sampe"
  else
    data[:mode] = "samse"
  end
  data[:tag] = ARGV.shift
  data[:inputs] = ["#{ARGV.shift}"]
  
  data[:inputs] << ["#{ARGV.shift}"] if "sampe" == data[:mode]
  
  groups << data
end

data = groups[index]

tag = data[:tag]

output = File.join(output_base,"#{index}.sam")

cmd = "bwa #{data[:mode]} -r \"#{tag}\" -f #{output} #{reference} #{data[:inputs].join(" ")}"

puts cmd
STDOUT.flush

exec cmd
