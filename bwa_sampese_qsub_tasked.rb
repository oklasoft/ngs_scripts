#!/usr/bin/env ruby1.8

args = ARGV.clone
# args.shift
output_base = args.shift
reference = args.shift
index = (ENV['SGE_TASK_ID']||1).to_i - 1

groups = [] #ARGV.each_slice(5).to_a[index]

while a = args.shift
  data = {}
  if "paired" == a then
    data[:mode] = "sampe"
  else
    data[:mode] = "samse"
  end
  data[:tag] = args.shift
  data[:inputs] = ["#{args.shift}","#{args.shift}"]
  
  data[:inputs] << ["#{args.shift}","#{args.shift}"] if "sampe" == data[:mode]
  
  groups << data
end

data = groups[index]

tag = data[:tag]

output = File.join(output_base,"#{index}.sam")

cmd = "bwa #{data[:mode]} -f #{output} #{reference} #{data[:inputs].join(" ")}"

puts cmd
STDOUT.flush

exec cmd
