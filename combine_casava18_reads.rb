#!/usr/bin/env ruby
#

require 'set'

base_reads = {} 

ARGV.each do |f|
  # file names will be like
  # SAMPLENAME_INDEX_L###_R#_###.fastq.gz
  sample_lane_read_name = File.basename(f,".fastq.gz").scan(/(.*_L\d+_R[1,2])_\d+/).first.first
  unless sample_lane_read_name
    raise "Can't figure out name for #{f}"
  end
  base_reads[sample_lane_read_name] ||= Set.new
  base_reads[sample_lane_read_name] << f
  puts "read name is: #{sample_lane_read_name}" if $DEBUG
  puts f if $DEBUG
end

base_reads.each do |base_lane_read,files|
  files = files.to_a.sort
  outdir = File.dirname(files.first)
  puts "#{base_lane_read} joined from #{files.join(", ")}"
  cmd = "zcat #{files.join(" ")} | gzip > #{outdir}/#{base_lane_read}.fastq.gz"
  system("qsub -cwd -V -m e -b y -j y -N #{base_lane_read} \"#{cmd}\"")
  puts cmd
end
