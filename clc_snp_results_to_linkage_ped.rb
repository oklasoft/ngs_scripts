#!/usr/bin/env ruby1.9

require 'csv'

GENOTYPE_IDS = {"A" => 1, "C" => 2, "G" => 3, "T" => 4}

unless ARGV.size >= 2
  STDERR.puts "Missing ped file & snp results file"
  exit 1
end

ped_file = ARGV.shift
snp_results_file = ARGV.shift

range_low = ARGV.shift
range_high = ARGV.shift

ped_data = {}
File.foreach(ped_file) do |line|
  parts = line.chomp.split(/\W+/)
  ped_data[parts[1]] = {:ped_id => parts[0], :subject => parts[1], :father => parts[2], :mother => parts[3], :sex => parts[4], :status => parts[5], :alleles => []}
end

STDERR.puts "Read in #{ped_data.size} subjects"

annotation_headers = ["contig", "position", "number of subjects", "reference", "allele variations"]
snps = []

multi_variation_snps = []

range_low = range_low.to_i if range_low
range_high = range_high.to_i if range_high

CSV.foreach(snp_results_file, {:headers => true}) do |row|
  subject_headers = row.headers - annotation_headers
  snp = "#{row.field('contig').gsub(/ /,'_')}_#{row.field('position').gsub(/,/,"")}"
  
  position = row.field('position').gsub(/,/,'').to_i
  
  next if range_low && position < range_low
  next if range_high && position > range_high
  
  variants = row.field('allele variations').split(/\/|,/).uniq
  if variants.size > 2  then
    multi_variation_snps << {:snp => snp, :variants => row.field('allele variations')}
    next
  end
  snps << snp
  alleles_seen = []
  
  subject_headers.each do |subject|
    # puts "adding snp #{snp} for #{subject}"
    ped_data[subject][:alleles] << row.field(subject)
    alleles_seen += row.field(subject).each_char.to_a
  end
  
  if alleles_seen.uniq.size > 2 then
    multi_variation_snps << {:snp => snp, :variants => alleles_seen.uniq.join("")}
    snps.pop
    subject_headers.each do |subject|
      ped_data[subject][:alleles].pop
    end
  end
  
end

STDERR.puts "#{snps.size} good snps found, #{multi_variation_snps.size} multi variant SNPs skipped"
STDERR.puts ""
STDERR.puts "#{snps.join(",")}"
STDERR.puts ""

tri_alleles = []

ped_data.each do |subject,data|
  print "#{data[:ped_id]} #{data[:subject]} #{data[:father]} #{data[:mother]} #{data[:sex]} #{data[:status]} "
  genotypes = []
  
  if data[:alleles].empty? then # a subject from the ped file for relationships
    genotypes = snps.map {|s| "0 0"}
  else
    data[:alleles].each_with_index do |g,i|
      # alleles = g.split(//)
      # STDERR.puts "allele: #{alleles.join(',')}"
      if 2 != g.length then # triallelic bastard
        genotypes << "0 0"
        tri_alleles << {:subject => subject, :snp => snps[i], :alleles => g}
      else
        genotypes << (g.each_char.map {|a| GENOTYPE_IDS[a]|| "0"}).join(" ")
      end
    end
  end
  raise "Error for #{subject} only #{genotypes.size}" if genotypes.size != snps.size
  puts "#{genotypes.join(" ")}"
end

unless multi_variation_snps.empty? then
  multi_variation_snps.each do |t|
    STDERR.puts t.inspect
  end
end

unless tri_alleles.empty? then
  tri_alleles.each do |t|
    STDERR.puts t.inspect
  end
end