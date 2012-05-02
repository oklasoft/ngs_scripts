#!/usr/bin/env ruby
#
# vaast_exon_extractor.rb
# Stuart Glenn 20120502
# Created by Stuart Glenn on 2011-06-07
#
# == Synopsis
# Quick script to extract some exon info from vaast output
#
# ==Author
# Stuart Glenn <Stuart-Glenn@omrf.org>
#
# ==Copyright
#  Copyright (c) 2011 Stuart Glenn, Oklahoma Medical Research Foundation. (OMRF)
#  All rights reserved.
#  Redistribution and use in source and binary forms, with or without
#  modification, are permitted provided that the following conditions are met:
#  1. Redistributions of source code must retain the above copyright
#     notice, this list of conditions and the following disclaimer.
#  2. Redistributions in binary form must reproduce the above copyright
#     notice, this list of conditions and the following disclaimer in the
#     documentation and/or other materials provided with the distribution.
#  3. All advertising materials mentioning features or use of this software
#     must display the following acknowledgement:
#     This product includes software developed by the OMRF
#  4. Neither the name of the Oklahoma Medical Research Foundation nor the
#     names of its contributors may be used to endorse or promote products
#     derived from this software without specific prior written permission.
#
#  THIS SOFTWARE IS PROVIDED BY COPYRIGHT HOLDERS AND CONTRIBUTORS ''AS IS'' AND ANY
#  EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
#  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
#  DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
#  DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
#  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
#  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
#  ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
#  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
#  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

OUTPUT_TOKEN = "\t"
output = []
ARGF.each do |line|
  next if line =~ /^#/
  if line =~ /^>/
    output += line.chomp.sub(/^>/,'').split(/\s+/)
    next
  elsif 2 == output.size && line !~ /.*:/
    (chr,position_str) = line.chomp.scan(/^(\w*)\s+[-+]\s+(.*)$/).first
    unless chr && !chr.empty? && position_str && !position_str.empty?
      raise "Error parsing for positions of #{output.join(",")} near line #{$.} working with '#{line.chomp}'"
    end
    output << chr
    positions = position_str.split(/\s/)
    #STDERR.puts "DEBUG: #{positions.inspect} from #{position_str}"
    output << positions.first.split(/;/).first
    output << positions.last.split(/;/)[1]
    
    puts output.join(OUTPUT_TOKEN)
    output = []
  end
end

puts output.join(OUTPUT_TOKEN) unless output.empty?
