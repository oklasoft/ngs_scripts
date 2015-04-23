#!/usr/bin/env ruby

class Person
  attr_reader :gid, :iid, :fid, :pid, :mid, :sex, :pheno
  def initialize(*parts)
    (@gid, @fid, @iid, @pid, @mid, @sex, @pheno) = parts
  end

  def to_s
    self.join("\t")
  end

  def join(sep="\t")
    [@fid,@iid,@pid,@mid,@sex,@pheno].join(sep)
  end
end

if ARGV.size < 3
  STDERR.puts "usage: #{ARGV[0]} [XCNV_INPUT] [PEDIGREE_TEMPLATE] [OUTPUT_BASE]"
  exit 1
end

xcnv_file = ARGV.shift
ped_template = ARGV.shift
output_base = ARGV.shift

puts "Converting #{xcnv_file} using #{ped_template} to #{output_base}"

ped = []
IO.foreach(ped_template) do |line|
  parts = line.chomp.split(/\t/)
  ped << Person.new(*parts)
end

used_ped = []
File.open("#{output_base}.cnv","w") do |out|
  IO.popen(["xcnv_to_cnv",xcnv_file]).each do |cnv_line|
    parts = cnv_line.split(/\t/)
    if parts[0] !~ /^FID/
      subject = ped.detect {|p| p.gid == parts[0]}
      if subject
        used_ped << subject
        parts[0] = subject.fid
        parts[1] = subject.iid
      end
    end
    out.puts parts.join("\t")
  end
end

File.open("#{output_base}.fam","w") do |out|
  used_ped.uniq.each do |p|
    out.puts p
  end
end

cmd = %W(plink --cnv-make-map --cnv-list #{output_base}.cnv --out #{output_base})
system(*cmd)
File.unlink("#{output_base}.log")
