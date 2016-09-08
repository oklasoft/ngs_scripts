#!/usr/bin/env ruby

$:.unshift File.dirname(File.realpath(__FILE__))
require 'tmpdir'
require 'tempfile'
require 'uri'
require 'uri-o3'
require 'optparse'
require 'shellwords'

MAX_DOWNLOAD_TRIES = 3

@options = {
  debug:false,
  verbose:false,
  tmp_base:Dir.tmpdir(),
  output_base:Dir.pwd(),
  reference:nil,
  trim:[],
  download_timeout:0,
  source_env:false,
  source_env_path:File.expand_path("~/.swiftrc"),
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

  o.on("--download-timeout INT",OptionParser::DecimalInteger,"Timeout file download after INT seconds, then retry, defaults to 0 (no timeout)") do |t|
    require 'timeout'
    @options[:download_timeout] = t
  end

  o.on("--source-env [FILE]","Source OS_AUTH_TOKEN from FILE shell env file, defaults to #{@options[:source_env_path]}") do |t|
    @options[:source_env] = true
    @options[:source_env_path] = File.expand_path(t) unless nil == t
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

env = {}
if @options[:source_env] && File.exist?(@options[:source_env_path])
  # This env parsing section extracted from https://github.com/bkeepers/dotenv
  # This module was MIT licensed Copyright (c) 2012 Brandon Keepers)
  LINE = /
        \A
        (?:export\s+)?    # optional export
        ([\w\.]+)         # key
        (?:\s*=\s*|:\s+?) # separator
        (                 # optional value begin
          '(?:\'|[^'])*'  #   single quoted value
          |               #   or
          "(?:\"|[^"])*"  #   double quoted value
          |               #   or
          [^#\n]+         #   unquoted value
        )?                # value end
        (?:\s*\#.*)?      # optional comment
        \z
      /x
  File.open(@options[:source_env_path],"rb:bom|utf-8").each do |line|
    line.chomp!
    if (match = line.match(LINE))
      key, value = match.captures
      next unless key == "OS_AUTH_TOKEN"
      env[key] = value.strip.sub(/\A(['"])(.*)\1\z/, '\2')
      break
    end
  end
end

delete_files = []
data[:inputs].map! do |i|
  if i =~ /^(https?|o3):\/\//
    u = URI(i)
    base = File.basename(u.path)
    f = Tempfile.new(base,@options[:tmp_base])
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
      env["OS_STORAGE_URL"] = u.os_storage_url
      %W/swift download -R 3 -o #{tmp_file} #{u.container} #{u.object}/
    end
    attempt = 0
    keep_trying = true
    while keep_trying && attempt < MAX_DOWNLOAD_TRIES do
      puts "Downloading (attempt #{attempt}) #{i}..."
      pid = spawn(env,*cmd,STDOUT=>STDERR,pgroup:true)
      status = nil
      begin
        Timeout::timeout(@options[:download_timeout]) do
        pid, status = Process.wait2(pid)
        keep_trying = false
        break
      end
      rescue Timeout::Error
        Process.kill(-15, pid)
        pid, status = Process.wait2(pid)
        attempt+=1
      end
    end
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

return_val = 0
Dir.mktmpdir(index.to_s,@options[:tmp_base]) do |tmp_prefix_dir|
  unless @options[:trim].empty?
    puts "Trimming" if @options[:verbose]
    if @options[:debug]
      @options[:trim].each do |t|
        puts "\t#{t}"
      end
    end
    trimmed_outputs = data[:inputs].map do |i|
      f = Tempfile.new("trimmo")
      f.close
      tmp_file = f.path + ".gz" #File.extname(base)
      f.unlink
      delete_files << tmp_file
      if "paired" == data[:mode]
        [tmp_file,"/dev/null"] #second is unpaired about which we do not care
      else
        tmp_file
      end
    end.flatten
    t = if "paired" == data[:mode]
          "trimmomatic_pe"
        else
          "trimmomatic_se"
        end
    cmd = "#{pre}#{t} -threads #{threads.to_i} -phred33 #{data[:inputs].join(" ")} #{trimmed_outputs.join(" ")} #{@options[:trim].join(" ")}#{post}"
    puts cmd
    STDOUT.flush
    unless @options[:debug]
      if system "/bin/bash", "-o", "pipefail", "-o", "errexit", "-c", cmd
        return_val = 0
      else
        return_val = $?.exitstatus
      end
    end
    if 0 != return_val
      $stderr.puts "Failed to trimmomatic"
    end
    data[:inputs] = trimmed_outputs.reject {|i| "/dev/null" == i}
    pre = ""
    post = ""
  end #if have trim
  if 0 == return_val
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
