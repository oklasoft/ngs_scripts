#!/usr/bin/env ruby

require 'optparse'

HUMAN_CHRS = %w/1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y MT/

MOUSE_CHRS = %w/1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 X Y M/.map{|c| "chr#{c}"}

chromosomes = HUMAN_CHRS

options = {:do_all => false, :memory => 16, :debug=> false, :extra_args=>[]}
optp = OptionParser.new
optp.banner = "Usage: #{File.basename(__FILE__)} "

optp.on("-r","--reference REFERENCE","Use the named REFERENCE file") do |conf|
  options[:reference] = File.expand_path(conf)
end

optp.on("-d","--dnsnp VCF","Use VCF file for filling rsID column") do |conf|
  options[:dbsnp] = File.expand_path(conf)
end

optp.on("-b","--output-base DIR","Save the output into DIR") do |conf|
  options[:output_base] = File.expand_path(conf)
end

optp.on("-p","--prefix PRE","Save the outout file using the PRE prefix") do |conf|
  options[:output_prefix] = conf
end

optp.on("-i","--input BAM","Process alignments from named input BAM file") do |conf|
  options[:input_bam] = File.expand_path(conf)
end

optp.on("-a","--all","Do all chromosomes at once instead of split by chr") do
  options[:do_all] = true
end

optp.on("-m","--memory GB", Integer,"Override default GATK memory setting of (#{options[:memory]}) in GB") do |conf|
  options[:memory] = conf.to_i
end

optp.on("-g","--genome ORGANISM",%w/human mouse/,"Select a genome set of chromosomes") do |org|
  if "mouse" == org
    chromosomes = MOUSE_CHRS
  end
end

optp.on("-E ARG", "Add extra argument ARG to GATK call") do |arg|
  options[:extra_args] << arg
end

optp.on("-D","--debug", "Debugging, print command but do not run") do
  options[:debug] = true
end

optp.on("-h","--help") do
  puts optp
  exit
end

optp.parse!
unless options[:output_base] && options[:output_prefix] && options[:input_bam] && options[:reference]
  STDERR.puts optp.help()
  exit(1)
end

ENV['JAVA_MEM_OPTS'] = "-Xmx#{options[:memory]}G"
threads = (ENV['NSLOTS']) || 1

cmd = %W/gatk -T HaplotypeCaller
         --pair_hmm_implementation VECTOR_LOGLESS_CACHING
         -ERC GVCF
         -nct #{threads}
         -R #{options[:reference]}
         -I #{options[:input_bam]}
         -A MQRankSum 
         -A ReadPosRankSum
         /
if options[:dbsnp] then
  cmd += %W/-D #{options[:dbsnp]}/
end

options[:extra_args].each do |ea|
  cmd += ea.split(/ /)
end

unless options[:do_all]
  # we will do it by chr
  begin
    Dir.mkdir(File.join(options[:output_base],'by_chr'))
  rescue
  end
  index = (ENV['SGE_TASK_ID']).to_i - 1
  chr = chromosomes[index]
  options[:output_prefix] = File.join(options[:output_base],"by_chr","#{options[:output_prefix]}-#{chr}")
  cmd += ["-L",chr]
else
  options[:output_prefix] = File.join(options[:output_base],"#{options[:output_prefix]}")
end
cmd += ["-o","#{options[:output_prefix]}.g.vcf.gz"]

puts cmd.join(" ")
STDOUT.flush

exec(*(cmd)) unless options[:debug]
