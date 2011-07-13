#!/usr/bin/env ruby1.9 -w -Ku
# encoding: UTF-8
#
# vcf_snp_cluster_summary.rb
# Created by Stuart Glenn on 2011-07-12T13:50:05-0500 
#
# == Synopsis
# Hack of script to read in the grepped out SNPCluster info & summarize those clusters
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
#

class VcfVariant
  attr_accessor :chromosome, :position, :rs_id, :reference, :observed_variant,
    :score, :tags
 
  def initialize(chr,pos,rs,ref,var,score,filters)
    @chromosome = chr
    @position = pos.to_i
    @rs_id = rs
    @reference = ref
    @observed_variant = var
    @score = score.to_f
    @tags = filters.split(/;/)
  end
  
  def distance_between(other)
    if @chromosome != other.chromosome
      return -1
    else
      return (@position - other.position).abs
    end
  end
  
  def gatk_standard?()
    @tags.include?("GATKStandard")
  end
  
end #VcfVariant class

class SnpCluster
  attr_reader :window_size
  def initialize(window)
    @window_size = window
    @snps = []
  end
  
  def empty?()
    @snps.empty?
  end
  
  def added_snp?(snp)
    if empty?()
      @snps << snp
      return true
    else
      @snps.each do |s|
        dist = s.distance_between(snp)
        return false if dist < 0 #a different chr
        if dist <= @window_size && dist >= 0
          @snps << snp
          return true
        end #if in cluster
      end #each snp
    end
    return false
  end #added_snp?
  
  def chromosome
    possible = @snps.map{|s| s.chromosome}.uniq
    raise "Too many chromosomes!" if possible.size > 1
    return possible.first
  end
  
  def size
    @snps.size
  end
  
  def start
    @snps.map{|s| s.position}.min
  end

  def stop
    @snps.map{|s| s.position}.max
  end
  
  def num_gatk_standard
    @snps.find_all{|s| s.gatk_standard?}.size
  end
  
  def dist_to_neighbors
    @snps.sort! {|a,b| a.position <=> b.position}
    dists = []
    @snps.each_with_index do |snp,i|
      dists << snp.distance_between(@snps[i-1]) unless i == 0
    end
    dists
  end
  
  # get some basic info like
  #how many SNPs make up the cluster, the start & stop of the cluster, 
  #how many of the SNPs are GATKStandard, 
  #and some basic summary info on distances between SNPs in that cluster
  def report(out)
    start = start()
    stop = stop()
    out.print "Num Snps:#{size()}"
    out.print "\tChromosome:#{chromosome()}"
    out.print "\tStart:#{start}"
    out.print "\tStop:#{stop}"
    out.print "\tTotal Size:#{stop-start}"
    out.print "\tGATK Standard:#{num_gatk_standard()}"
    out.print "\tDistances:#{dist_to_neighbors().join(";")}"
    out.puts
  end
  
end #SnpCluster class

clusters = []
current_cluster = SnpCluster.new(10)

ARGF.each do |line|
  next if line =~ /^#/
  snp = VcfVariant.new(*line.chomp.split(/\s+/))
  next if current_cluster.added_snp?(snp)
  
  current_cluster.report($stdout)
  
  clusters << current_cluster
  current_cluster = SnpCluster.new(10)
  current_cluster.added_snp?(snp)
end
current_cluster.report($stdout)

puts "Total number clusters: #{clusters.size}"