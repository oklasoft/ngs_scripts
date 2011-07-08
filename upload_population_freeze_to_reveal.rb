#!/usr/bin/env ruby1.9 -w -Ku
# encoding: UTF-8
#
# upload_population_freeze_to_reveal.rb
# Created by Stuart Glenn on 2011-07-07T15:30:10-0500
#
# == Synopsis
# Quicker hack to do work of "uploading" a freeze to the REVEAL system
# That is it will give you SQL to create interval records & copy versions
# of the input files
#
# == Inputs
#  - A file listing loci info of the target VCF(s): chr start stop name
#  - A reveal id of the population
#  - A reveal id of the freeze
#  - One or more directories named by loci of the for file types
#    - These file types are loci.map, loci.bed, loci.ped loci_summary.txt
#
# == Usage
#  upload_population_freeze_to_reveal.rb -l FILE -f ID -p ID -t DIR INPUT_DIR(s)
#
#  For help use upload_population_freeze_to_reveal.rb -h
#
# == Options
#  -h, --help             Display this help message
#  -V, --verbose          Increased verbosity of output
#  -l, --loci FILE        A file with information about possible loci
#  -f, --freeze id        A REVEAL id for the targeted freeze
#  -p, --populaation id   A REVEAL id for the targeted population
#
# ==Author
#  Stuart Glenn <Stuart-Glenn@omrf.org>
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

require 'optparse'
require 'ostruct'
require 'fileutils'

class RevealUploader
  
  CHR_IDS = {
    "chr1" => 1,
    "chr2" => 2,
    "chr3" => 3,
    "chr4" => 4,
    "chr5" => 5,
    "chr6" => 6,
    "chr7" => 7,
    "chr8" => 8,
    "chr9" => 9,
    "chr10" => 10,
    "chr11" => 11,
    "chr12" => 12,
    "chr13" => 13,
    "chr14" => 14,
    "chr15" => 15,
    "chr16" => 16,
    "chr17" => 17,
    "chr18" => 18,
    "chr19" => 19,
    "chr20" => 20,
    "chr21" => 21,
    "chr22" => 22,
    "chrX" => 23,
    "chrY" => 24,
  }
  
  EXTENSIONS = {
    :summary => "_summary.txt",
    :ped => ".ped",
    :map => ".map",
    :bed => ".bed"
  }
  
  def initialize(args,ios = {})
    @args = args
    set_inputs_outputs(ios)
    set_default_options()
  end
  
  def run
    unless options_parsed?
      @stderr.puts "Error parsing arguments"
      output_usage(@stderr)
      return nil
    end
    unless options_valid?
      @stderr.puts "Invalid options: #{@error_message}"
      output_usage(@stderr)
      return nil
    end

    had_errors = false

    @options.input_dirs.each do |dir_path|
      unless process_dir(dir_path)
        @stderr.puts "Error processing #{dir_path}: #{@error_message}"
        had_errors = true
      end
    end #each dir
    return !had_errors
  end #run
  
  
  private
  
  def process_dir(input_dir)
    locus = find_locus(File.basename(input_dir))
    unless locus
      @error_message = "Unable to determine locus name"
      return false
    end
    
    files = find_files_in(input_dir,locus[:name])
    unless found_right_files?(files)
      @error_message = "Unable to find all the right input files for #{locus[:name]}"
      return false
    end
    
    return add_locus_to_datastore(locus,files)
  end
  
  def add_locus_to_datastore(locus,files)
    store_files(files,locus[:name]) && insert_locus_to_db(locus)
  end
  
  def store_files(files,locus_name)
    dest_dir = File.join(@options.output_dir,locus_name)
    if Dir.exists?(dest_dir)
      @error_message = "The destination for #{locus_name} exists already at #{dest_dir}, skipping"
      return false
    end
    Dir.mkdir(dest_dir)
    files.each do |key,path|
      FileUtils.copy(path,File.join(dest_dir,"#{key}.txt"))
    end
    return true
  end
  
  def insert_locus_to_db(locus)
    sql = <<-EOF
INSERT INTO intervals(start,stop,locus,chromosome_id,race_id,release_id) VALUES(#{locus[:start]},#{locus[:stop]},'#{locus[:name]}',#{CHR_IDS[locus[:chr]]},#{@options.population},#{@options.freeze_id});
EOF
    $stdout.puts sql
    return true
  end
  
  def find_files_in(dir,locus_name)
    files = nil
    Dir.chdir(dir) do
      files = {}
      files[:summary] = Dir.glob("*-#{locus_name}_summary.txt")
      return nil unless files[:summary] && 1 == files[:summary].size
      files[:summary] = files[:summary].first
      prefix = File.basename(files[:summary],"_summary.txt")

      EXTENSIONS.each do |key,ext|
        files[key] = Dir.glob("#{prefix}#{ext}")
        return nil unless files[key] && 1 == files[key].size
        files[key] = File.join(dir,files[key].first)
      end
    end
    files
  end
  
  def found_right_files?(files)
    return false unless files && files.class == Hash
    return false unless files[:bed] && files[:map] && files[:ped] && files[:summary]
    return true
  end
  
  def find_locus(filename)
    @loci.each do |locus_name,locus_data|
      if filename =~ /^#{locus_name}$/
        return locus_data
      end
    end
    return nil
  end #find_locus
  
  def load_loci()
    if nil == @options.loci_file
      @error_message = "Missing loci file option"
      return false
    end
    @loci = {}
    IO.foreach(@options.loci_file) do |locus_line|
      next if locus_line =~ /^#/ || locus_line =~ /^$/
      (chr,start,stop,name) = locus_line.chomp.split(/\s+/)
      @loci[name] = {:chr => chr, :start => start, :stop => stop, :name => name}
    end
    return true
  rescue  => error
    @error_message = "#{error.message}"
    return false
  end #load_loci
  
  def set_inputs_outputs(ios)
    @stdin = ios[:stdin] || $stdin
    @stdout = ios[:stdout] || $stdout
    @stderr = ios[:stderr] || $stderr
  end

  def options_valid?
    load_loci() &&
    output_dir_valid?() &&
    freeze_given?() &&
    population_given?() &&
    input_dirs_given?()
  end #options_valid?
  
  def input_dirs_given?()
    unless @options.input_dirs
      @error_message = "Missing input dirs to process"
      return false
    end
    trouble = []
    @options.input_dirs.each_with_index do |dir,i|
      dir = File.expand_path(dir)
      @options.input_dirs[i] = dir
      unless File.directory?(dir)
        trouble << "Specified input dir (#{dir}) is not a directory"
        next
      end
      unless File.readable?(dir)
        trouble << "Specified input dir (#{dir}) is not readable"
      end
    end #each dir
    unless trouble.empty?
      @error_message = "Trouble with input dir(s): #{trouble.join("; ")}"
      return false
    end
    return true
  end
  
  def freeze_given?()
    id_given?(@options.freeze_id,"freeze")
  end

  def population_given?()
    id_given?(@options.population,"population")
  end
  
  def id_given?(test_id,name)
    if nil == test_id then
      @error_message = "Must specify the #{name}"
      return false
    end
    unless test_id.to_i > 0 then
      @error_message = "Unknown #{name} of #{test_id} given"
      return false
    end
    return true
  end
  
  def output_dir_valid?
    unless @options.output_dir
      @error_message = "Must specify a base output directory"
      return false
    end
    @options.output_dir = File.expand_path(@options.output_dir)
    unless File.directory?(@options.output_dir)
      @error_message = "Given (#{@options.output_dir}) output directory is NOT a directory"
      return false
    end
    unless File.writable?(@options.output_dir)
      @error_message = "Given (#{@options.output_dir}) output directory is not writable"
      return false
    end
    return true
  end

  def set_default_options()
    @options = OpenStruct.new(
      :verbose => false,
      :output_dir => nil,
      :freeze_id => nil,
      :population => nil,
      :loci_file => nil,
      :input_dirs => nil
    )
  end #set_default_options

  def options_parsed?
    opt_parser = OptionParser.new() do |opts|
      opts.on('-h','--help') { output_help(@stdout); exit(0) }
      opts.on('-V', '--verbose')    { @options.verbose = true }

      opts.on("-l","--loci", "=REQUIRED") do |file|
        @options.loci_file = file
      end

      opts.on("-o","--output", "=REQUIRED") do |file|
        @options.output_dir = file
      end

      opts.on("-f","--freeze", "=REQUIRED") do |id|
        @options.freeze_id = id
      end

      opts.on("-p","--population", "=REQUIRED") do |id|
        @options.population = id
      end
    end

    opt_parser.parse!(@args) rescue return false
    unless @args.empty?
      @options.input_dirs = @args.find_all {|a| File.directory?(a)} 
    end
    return true
  end #options_parsed?
  
  def output_usage(out)
    out.puts <<-EOF
#{File.basename($0)} -l FILE -f ID -p ID -o OUTPUT_BASE_DIR INPUT_DIR(s)

-h, --help             Display this help message
-V, --verbose          Increased verbosity of output
-l, --loci FILE        A file with information about possible loci
-f, --freeze id        A REVEAL id for the targeted freeze
-p, --population id   A REVEAL id for the targeted population
-o, --output DIR       A base output target dir for the new files
EOF
  end
  
end

if $0 == __FILE__
  stat = RevealUploader.new(ARGV.clone).run
  if true == stat
    exit 0
  else
    $stderr.puts "Errors in running, check output" unless nil == stat
    exit 1
  end
end