#!/usr/bin/env ruby
#
# Stuart Glenn <Stuart-Glenn@omrf.org>
# 2011-02-07
# 
# Copyright (c) 2011, Oklahoma Medical Research Foundation
# All rights reserved.
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#     * Redistributions in binary form must reproduce the above copyright
#       notice, this list of conditions and the following disclaimer in the
#       documentation and/or other materials provided with the distribution.
#     * Neither the name of the Oklahoma Medical Research Foundation nor the
#       names of its contributors may be used to endorse or promote products
#       derived from this software without specific prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL OKLAHOMA MEDICAL RESEARCH FOUNDATION BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# combine GATK DepthOfCoverage outputs

BASES = %w/sample_cumulative_coverage_counts
sample_cumulative_coverage_proportions
sample_interval_statistics
sample_interval_summary
sample_statistics
sample_summary/

if ARGV.size < (BASES.size + 2)
  STDERR.puts "Usage: #{$0} <output_base> <input_files>"
  exit 1
end

output_base = ARGV.shift
inputs = ARGV.clone

samples = {}

inputs.each do |input|
  file = File.basename(input)
  parts = file.split(/\./)
  id = parts[0]
  if 1 == parts.size
    raise "Duplicate sample id #{id}" if samples[id]
    samples[id] = {:full => input}
  else
    next unless BASES.include?(parts[1])
    samples[id][parts[1].to_sym] = input
  end
end

# sample_summary
# get one header line, then add each line, final line is Total, sum on col 2, sum col 3, then N/A
File.open("#{output_base}.sample_summary","w") do |out|
  print "sample_summary: "
  header_done = false
  total_total = 0.0
  total_avg = 0.0
  
  samples.each do |sample,files|
    print "."
    lines = IO.readlines(files[:sample_summary])    
    unless header_done
      out.puts lines[0].chomp
      header_done = true
    end
    out.puts lines[1].chomp
    parts = lines[1].chomp.split(/\t/)
    total_total += parts[1].to_f
    total_avg += parts[2].to_f
  end
  out.puts ["Total",total_total,total_avg,"N/A","N/A","N/A"].join("\t")
  puts ""
end

# sample_statistics
# just copy the lines, its an easy one
File.open("#{output_base}.sample_statistics","w") do |out|
  print "sample_statistics: "
  header_done = false
  
  samples.each do |sample,files|
    print "."
    lines = IO.readlines(files[:sample_statistics])    
    unless header_done
      out.puts lines[0].chomp
      header_done = true
    end
    out.puts lines[1].chomp
  end
  puts ""
end

# sample_interval_summary
# target, sum cols 1, sum cols2, join column sets of 6 at a time
File.open("#{output_base}.sample_interval_summary","w") do |out|
  print "sample_interval_summary: "
  header = []
  targets = []
  total_coverages = []
  total_averages = []
  target_lines = []
  
  samples.keys.each_with_index do |sample,i|
    print "."
    IO.foreach(samples[sample][:sample_interval_summary]) do |line|
      parts = line.chomp.split(/\t/)
      if 1 == $. then
        # header, if first line all, if not first toss first two columns
        if 0 != i then
          parts.shift
          parts.shift
          parts.shift
        end
        header += parts
      else
        # data lines
        target = parts.shift
        coverage = parts.shift.to_i
        average = parts.shift.to_f
        if 0 == i then
          targets << target
          total_coverages << coverage
          total_averages << average
          target_lines << parts
        else
          index = $. - 2
          total_coverages[index] += coverage
          total_averages[index] += average
          target_lines[index] += parts
        end  # if first
      end # if header
    end #each line
  end #each sample
  puts ""

  out.puts header.join("\t")
  targets.each_with_index do |t,i|
    out.puts "#{t}\t#{total_coverages[i]}\t#{total_averages[i]}\t#{target_lines[i].join("\t")}"
  end  
end


# sample_interval_statistics
# copy header each col shift max
File.open("#{output_base}.sample_interval_statistics","w") do |out|
  print "sample_interval_statistics: "
  depths = []
  samples.keys.each_with_index do |sample,i|
    print "."
    lines = IO.readlines(samples[sample][:sample_interval_statistics])
    if 0 == i then
      out.puts lines[0].chomp
      (lines[1].split(/\t/).size-1).times {depths << []}
    end
    parts = lines[1].split(/\t/)
    parts.shift
    parts.each_with_index do |p,i|
      depths[i] << p.to_i
    end
  end #each sample
  
  depths.each do |d|
    d.sort!
    d.reverse!
  end
  
  samples.keys.each_with_index do |sample,i|
    out.print "At_least_#{i+1}_samples"
    depths.each do |d|
      out.print "\t#{d.shift}"
    end
    out.puts
  end
  
  puts ""
end


# sample_cumulative_coverage_proportions
# just copy the lines, its an easy one
File.open("#{output_base}.sample_cumulative_coverage_proportions","w") do |out|
  print "sample_cumulative_coverage_proportions: "
  header_done = false
  
  samples.each do |sample,files|
    print "."
    lines = IO.readlines(files[:sample_cumulative_coverage_proportions])    
    unless header_done
      out.puts lines[0].chomp
      header_done = true
    end
    out.puts lines[1].chomp
  end
  puts ""
end

# sample_cumulative_coverage_counts
# copy header each col shift max?
File.open("#{output_base}.sample_cumulative_coverage_counts","w") do |out|
  print "sample_cumulative_coverage_counts: "
  depths = []
  samples.keys.each_with_index do |sample,i|
    print "."
    lines = IO.readlines(samples[sample][:sample_cumulative_coverage_counts])
    if 0 == i then
      out.puts lines[0].chomp
      (lines[1].split(/\t/).size-1).times {depths << []}
    end
    parts = lines[1].split(/\t/)
    parts.shift
    parts.each_with_index do |p,i|
      depths[i] << p.to_i
    end
  end #each sample
  
  depths.each do |d|
    d.sort!
    d.reverse!
  end
  
  samples.keys.each_with_index do |sample,i|
    out.print "NSamples_#{i+1}"
    depths.each do |d|
      out.print "\t#{d.shift}"
    end
    out.puts
  end
  
  puts ""
end

# 'all'
# target, sum of row, average of row, join column sets of 1
File.open("#{output_base}","w") do |out|
  print "full: "
  header = []
  targets = []
  total_depth = []
  target_lines = []
  
  keys = samples.keys
  positions = Array.new(keys.size - 1,0)
  header_done = false

  IO.foreach(samples[keys.shift][:full]) do |first_file_line|
    print "." if (0 == $. % 100000)
    new_line = first_file_line.chomp.split(/\t/) # we always start the first columns of first
    new_line[1] = new_line[1].to_i if header_done
    keys.each_with_index do |k,i|
      File.open(samples[k][:full]) do |fin|
        fin.seek(positions[i],IO::SEEK_SET)
        (target,depth,avg,data) = fin.readline.chomp.split(/\t/)
        positions[i] = fin.pos
        depth = depth.to_i
        new_line[1] += depth if header_done
        new_line << data
      end
    end
    
    new_line[2] = new_line[1].to_f/samples.keys.size if header_done
    out.puts new_line.join("\t")
    header_done = true
  end

end

