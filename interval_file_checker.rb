#!/usr/bin/env ruby


class Interval
  attr_accessor :chr, :start, :stop
  def initialize(chr,start,stop)
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
  merged_intervals = []
  ARGF.each do |line|
    (chr,start,stop) = line.chomp.scan(/(.*):(\d+)-(\d+)/).first
    next unless chr && start && stop
    i = Interval.new(chr,start,stop)
    merged_intervals.each do |mi|
      if mi.overlaps?(i)
        mi.merge!(i)
        i = nil
      end
    end
    merged_intervals << i if i
  end

  merged_intervals.sort! do |a,b|
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
  puts merged_intervals.join("\n")
#end
