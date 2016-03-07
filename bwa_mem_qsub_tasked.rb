#!/usr/bin/env ruby

require 'tmpdir'
require 'tempfile'
require 'uri'
require 'optparse'

module URI
  class O3 < Generic
    USE_REGISTRY = true
  end
  @@schemes['O3'] = O3
end

@options = {
  debug:false,
  verbose:false,
  tmp_base:Dir.tmpdir(),
  output_base:Dir.pwd(),
  reference:nil,
  trim:[],
}

op = OptionParser.new do |o|
  o.banner = "Usage: #{File.basename(__FILE__)} -r REF [-t TMP_BASE] [-o OUTPUT_BASE] "

  o.on("-r","--reference FILE","Use FILE as alignment reference to bwa") do |r|
    @options[:reference] = File.expand_path(r)
  end
  o.on("-t","--tmp DIR","Use DIR as base directory for temp working files, defaults to #{@options[:tmp_base]}") do |t|
    @options[:tmp_base] = File.expand_path(t)
  end
  o.on("-o","--out DIR","Use DIR as base directory for final output, defaults to cwd") do |t|
    @options[:output_base] = File.expand_path(t)
  end

  o.on("--trim STRING","Trim fastq first using trimmomatic with STRING, can be specified multiple times") do |t|
    @options[:trim] << t
  end

  o.on("-v","--verbose","Increase verbosity of output") do
    @options[:verbose] = true
  end
  o.on("-D","--debug","Enable debugging mode, does not actually upload") do
    @options[:debug] = true
  end
  o.on("-h","--help","Show this help message") do
    puts o
    exit(0)
  end
  o.parse!
end

unless @options[:reference]
  $stderr.puts "Missing alignment reference"
  $stderr.puts op.help()
  exit(1)
end

args = ARGV.clone
index = (ENV['SGE_TASK_ID']||1).to_i - 1
threads = (ENV['NSLOTS']) || 1

groups = []

while a = args.shift
  data = {}
  if "paired" == a then
    data[:mode] = "paired"
  else
    data[:mode] = "single"
  end
  data[:tag] = args.shift
  data[:inputs] = ["#{args.shift}"]
  data[:inputs] << "#{args.shift}" if "paired" == data[:mode]

  groups << data
end

data = groups[index]
tag = data[:tag]
output = File.join(@options[:output_base],"#{index}.bam")

delete_files = []
data[:inputs].map! do |i|
  if i =~ /^(https?|o3):\/\//
    u = URI(i)
    base = File.basename(u.path)
    f = Tempfile.new(base)
    f.close
    tmp_file = f.path + File.extname(base)
    f.unlink
    delete_files << tmp_file
    cmd = case u.scheme
    when /^https?/
      # curl it
      %W/curl -s --retry 3 -f -o #{tmp_file} #{i}/
    when /^o3/
      # swift it
      # TODO pull out account with .registry & create the OS_STORAGE_URL env
      (container,object)= u.path.scan(/^\/([^\/]+)\/(.*)/)[0]
      %W/swift download -R 3 -o #{tmp_file} #{container} #{object}/
    end
    puts "Downloading #{i}..."
    pid = spawn(*cmd,STDOUT=>STDERR)
    pid, status = Process.wait2(pid)
    if nil == status
      raise "Unable to start object download"
    elsif 0 != status.exitstatus
      raise "Failure downloading object: #{$?}"
    end
    tmp_file
  else
    i
  end
end

pre = ""
post = ""
data[:inputs].map! do |i|
  if i =~ /\.xz$/
    pre ="(t=`mktemp`;"
    post = " && rm ${t})"
    "<(xzcat #{i} || rm -f ${t})"
  else
    i
  end
end

return_val = -1
Dir.mktmpdir(index.to_s,@options[:tmp_base]) do |tmp_prefix_dir|
  unless @options[:trim].empty?
    puts "Trimming" if @options[:verbose]
    if @options[:debug]
      @options[:trim].each do |t|
        puts "\t#{t}"
      end
    end
  end
  cmd = "#{pre}bwa mem -v 1 -M -t #{threads.to_i-2} -R \"#{tag}\" #{@options[:reference]} #{data[:inputs].join(" ")}#{post}"
  cmd += "|samtools view -Shu - | samtools sort -O bam -@ 4 -m 4G -T #{tmp_prefix_dir}/#{index} > #{output}"

  puts cmd
  STDOUT.flush
  unless @options[:debug]
    if system "/bin/bash", "-o", "pipefail", "-o", "errexit", "-c", cmd
      return_val = 0
    else
      return_val = $?.exitstatus
    end
  end
end

if 0 == return_val
  delete_files.each do |f|
    begin
      File.delete(f) if File.exists?(f)
    rescue
    end
  end
end
exit return_val
