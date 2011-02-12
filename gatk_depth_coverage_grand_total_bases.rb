#!/usr/bin/env ruby
# 
# Stuart Glenn <Stuart-Glenn@omrf.org> 20110212
#
# #
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

TOTAL_INDEX = 1
FIRST_SAMPLE_INDEX = 3

total = 0
totals = Hash.new(0)
samples = []

input = nil
if ARGV.size > 0 then
  input = File.open(ARGV.shift)
else
  input = $stdin
end

input.each_line do |line|
  parts = line.chomp.split(/\t/)
  
  if 1 == $.
    parts[FIRST_SAMPLE_INDEX,parts.size].each do |sample|
      sample.sub!(/Depth_for_/,'')
      samples << sample.to_sym
    end
  else
    total += parts[TOTAL_INDEX].to_i
    samples.each_with_index do |sample,i|
      totals[sample] += parts[i+FIRST_SAMPLE_INDEX].to_i
    end
  end
end

input.close

puts "Sample\tTotal Bases"
samples.each do |sample|
  puts "#{sample}\t#{totals[sample]}"
end
puts "TOTAL\t#{total}"
