#!/usr/bin/env ruby


class Interval
  attr_accessor :chr, :start, :stop
  def initialize(chr,start,stop)
    raise ArgumentError, "Missing values" unless chr && start && stop && !start.empty? && !chr.empty? && !stop.empty?
    @chr = case chr
           when "X"
             23
           when "Y"
             24
           else
             chr.to_i
           end
    @start = start.to_i
    @stop = stop.to_i
  end

  def size
    @stop - @start
  end

  def to_s
    chr = case @chr
          when 23
            "X"
          when 24
            "Y"
          else
            @chr
          end
    "#{chr}:#{@start}-#{@stop}"
  end

  def overlaps?(b)
    return false unless self.chr == b.chr
    return true if self.start <= b.start && self.stop >= b.start
    return true if self.start <= b.stop && self.stop >= b.stop
    return true if b.start <= self.stop && b.stop >= self.stop
    return true if b.start <= self.start && b.stop >= self.start
  end

  def merge!(b)
    raise ArgumentError.new("Must be same chromosome") unless self.chr == b.chr
    self.start = [self.start,b.start].min
    self.stop = [self.stop,b.stop].max
  end
end

#if __FILE__ == ARGV[0]
  merged_intervals = Hash.new {|hash,key| hash[key] = []}
  ARGF.each do |line|
    line.chomp!
    (chr,start,stop) = line.scan(/(.*):(\d+)-(\d+)/).first
    i = Interval.new(chr.chomp,start.chomp,stop.chomp)
    if nil == i
      raise "Failed to make interval from '#{line}' had chr=#{chr} start=#{start} stop=#{stop}"
    end
    merged_intervals[i.chr].each do |mi|
      if mi.overlaps?(i)
        mi.merge!(i)
        i = nil
        break
      end
    end
    merged_intervals[i.chr] << i if i
  end

  merged_intervals.keys.sort.each do |chr|
    merged_intervals[chr].sort! do |a,b|
      if a.chr == b.chr
        if a.start == b.start
          a.stop <=> b.stop
        else
          a.start <=> b.start
        end
      else
        a.chr.to_i <=> b.chr.to_i
      end
    end
    puts merged_intervals[chr].join("\n")
  end
#end
