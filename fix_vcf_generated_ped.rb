#!/usr/bin/env ruby1.9
# encoding: UTF-8
#
# fix_vcf_generated_ped.rb
# Created by Stuart Glenn on 2011-01-04
#
# == Synopsis
# Quick script to take a template ped file (first 6 columns) of real data &
# replace then then the bogus columns in vcftools produced ped files. We will
# also optionally add invidiuals to the new ped files from the template. These
# added individuals will have no genotypes for all the genotypes listed in the
# original
#
# == Inputs
#  - The template file is 7 columns, first being indivual id in the target ped files
#  followed then by the standard 6 ped columns of plink, all separated by tabs
#  family individual father mother sex(0=unknown,1=male,2=female) phenotype(-9,0=missing,1=unaffected,2=affected)
#  - The rest of the files will be plain text ped files with the first 6 columns being 
#  ped info, then followed by any genotypes
#
# == Usage
#  fix_vcf_generated_ped.rb -t PED_TEMPLATE -p OUTPUT_PREFIX INPUT_PEDIGREE(S) 
#
#  For help use fix_vcf_generated_ped.rb -h
#
# == Options
#  -h, --help             Display this help message
#  -v, --version          Display the version information
#  -V, --verbose          Increased verbosity of output
#  -t, --template FILE    Specify the input template pedigree
#  -a, --add_indivuals    Add any/all individuals in the template, regardless of them being present in original
#  -p, --prefix STR       Specify the prefix the output file(s)
#  -i, --inplace EXT      Save new file in place, saving original files with given extension.
#                         Note original is not saved if no extension is given with -i option
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
require 'tempfile'
require 'fileutils'

class PedFixerApp
  VERSION       = "0.0.9-pre01"
  REVISION_DATE = "2011-01-04"
  AUTHOR        = "Stuart Glenn <Stuart-Glenn@omrf.org>"
  COPYRIGHT     = "Copyright (c) 2011 Oklahoma Medical Research Foundation"
  
  NUM_PEDIGREE_TEMPLATE_FIELDS = 7
  
  NUM_PEDIGREE_FIELDS = 6
  
  PEDIGREE_FILE_DELIMITER = "\t"
  
  # Create an app object to run
  # * +args+ - Array of command line argument string
  # * +ios+ - Optional hash of :stdin, :stdout, :stderr to set to corresponding @std
  def initialize(args,ios = {})
    @args = args
    set_inputs_outputs(ios)
    set_default_options()
  end

  # Run the application
  def run
    if options_parsed? && options_valid?
      unless load_pedigree_template()
        @stderr.puts("")
        output_usage(@stderr)
        exit(1)
      end

      @options.input_files.each do |input|
        fix_ped_file(input)
      end
    else
      @stderr.puts("")
      output_usage(@stderr)
      exit(1)
    end
  end #run
  
  
  private
  
  # Read the ped template file & load its data into a hash for later use
  def load_pedigree_template
    @ped_template = {}
    @individuals_not_seen = []
    IO.foreach(@options.template_file) do |line|
      next if line =~ /^#/
      parts = line.chomp.split(PEDIGREE_FILE_DELIMITER)
      unless NUM_PEDIGREE_TEMPLATE_FIELDS == parts.size
        @stderr.puts "Failure near line #{$.} of #{@options.template_file}; incorrect number of fields"
        return false
      end
      id = parts.shift
      @ped_template[id] = parts
      @individuals_not_seen << id
    end
    return true
  end #load_pedigree_template

  # Actually get to the work of fixing & saving a new file
  # * +input_file+ - The path to the input file to fix (or @stdin)
  # We will save the results based on what output_for_input() tells us
  def fix_ped_file(input_file)
    output = output_for_input(input_file)
    unless output_is_valid?(output)
      @stderr.puts "Invalid output for '#{input_file.inspect}' as #{output.inspect}"
      return false
    end
    
    output = File.new(output,"w") if String == output.class
    write_fixed_file_to(input_file,output)
    output.close
    
    finalize_output(input_file,output)
    return true
  end #fix_ped_file(input_file)
  
  # Does any last minute work on the the output, which mainly means fixing
  # the inplace naming changing
  # * +input+ - the original input file path
  # * +output+ - the output File object
  def finalize_output(input,output)
    return if @stdout == output || !@options.save_inplace
    if @options.backup_extension then
      swap_input_and_output_files(input,output.path)
    else
      move_temp_file_to_input!(output,input)
    end
  end #finalize_output(output)
  
  # Given two file paths, swap them
  def swap_input_and_output_files(src_path,dest_path)
    tmp_path = tmp_rename_file(src_path)
    File.rename(src_path,tmp_path)
    File.rename(dest_path,src_path)
    File.rename(tmp_path,dest_path)
  end #swap_input_and_output_files()
  
  # Get a short temp file for the rename
  def tmp_rename_file(input)
    f = Tempfile.new("tmp",File.dirname(input))
    path = f.path
    f.close
    f = nil
    return path
  end #tmp_rename_file(input)
  
  # Really copy the temp file over the input file
  # * This will destroy the original input
  # * +source_file+ - The original output file object, to get the source path
  # * +dest_path+ - What was the original input path, is now the destination
  def move_temp_file_to_input!(source_file,dest_path)
    FileUtils.copy(source_file.path,dest_path)
  end #move_temp_file_to_input(output,input)
  
  # Reads through the input file, fixing each line & writing it to the output
  # * +input+ - The path to the input pedigree
  # * +output+ - An opened & writable File object into which the results go
  # We will also use the @ped_template hash to try to fix each line
  def write_fixed_file_to(input,output)
    num_markers = 0
    IO.foreach(input) do |line|
      parts = line.chomp.split(PEDIGREE_FILE_DELIMITER)
      output.puts(fixed_line!(parts))
      if 0 == num_markers
        num_markers = (parts.size)/2
      end
    end
    add_missing_individuals(output,num_markers) if @options.add_individuals
  end #write_fixed_file_to(input,output)

  # Add any subjects from the template which we didn't yet add to output by 
  # giving them missing genotypes for everything
  # * +output+ - To where to write
  # * +num_markers+ - How many genotypes to are they to get
  def add_missing_individuals(output,num_markers)
    genos = ["0"] * num_markers * 2
    @individuals_not_seen.each do |id|
      ped = @ped_template[id]
      output.puts (ped + genos).join(PEDIGREE_FILE_DELIMITER)
    end
  end #add_missing_individuals(output)

  # Fix the given line by looking for the match in the ped_template hash
  def fixed_line!(line_parts)
    orig_ped = line_parts.shift(NUM_PEDIGREE_FIELDS)
    ped = @ped_template[orig_ped[0]] || orig_ped
    @individuals_not_seen.delete(orig_ped[0])
    (ped + line_parts).join(PEDIGREE_FILE_DELIMITER)
  end #fixed_line(line)

  # Setup the streams to use for stdin, out & error
  # * +ios+ - Hash with keys of :stdin, :stdout, :stderr
  # Any key not present the correspodning stream will be set the standard
  def set_inputs_outputs(ios)
    @stdin = ios[:stdin] || STDIN
    @stdout = ios[:stdout] || STDOUT
    @stderr = ios[:stderr] || STDERR
  end
  
  # Just make sure the app is happy with some steady & known default options
  # build the @options hash
  def set_default_options()
    @options = OpenStruct.new(
      :template_file => nil,
      :prefix => nil,
      :verbose => false,
      :save_inplace => false,
      :backup_extension => nil,
      :add_individuals => false,
      :input_files => [@stdin]
    )
  end #set_default_options
  
  # Do the work of parsing out the options in the @arg array into the right
  # places into the @options
  # ==== Retuns
  # * +true+ - For good option parsing
  # * +false+ - On errors during parsing
  def options_parsed?
    opts = OptionParser.new() do |opts|
      opts.on('-v','--version') { output_version(@stdout); exit(0) }
      opts.on('-h','--help') { output_help(@stdout); exit(0) }
      opts.on('-V', '--verbose')    { @options.verbose = true }

      opts.on('-a', '--add_individuals')    { @options.add_individuals = true }
      
      opts.on("-t","--template", "=REQUIRED") do |template_file|
        @options.template_file = template_file
      end

      opts.on("-p","--prefix", "=REQUIRED") do |prefix|
        @options.prefix = prefix
      end
      
      opts.on("-i","--inplace [EXTENSION]") do |ext|
        @options.save_inplace  = true
        @options.backup_extension = ext || nil
        @options.backup_extension.sub!(/\A\.?(?=.)/, ".") if @options.backup_extension
      end
    end
    
    opts.parse!(@args) #rescue return false
    @options.input_files = @args unless @args.empty?
    return true
  end #options_parsed?
  
  # Make sure the options we have are correct & valid
  def options_valid?
    template_valid?() &&
    output_options_valid?()
  end #options_valid?
  
  # Test the input template by checking it is readable & contains the valid 7 
  # columns
  def template_valid?()
    msg = ''
    if @options.template_file && File.readable?(@options.template_file)
      msg = "Wrong format, no data?"
      IO.foreach(@options.template_file) do |line|
        next if line =~ /^#/
        if NUM_PEDIGREE_TEMPLATE_FIELDS == line.split(PEDIGREE_FILE_DELIMITER).size
          return true
        else
          msg = "Wrong number of fields"
        end
      end
    else
      msg = "Unable to read file"
    end
    @stderr.puts "Pedigree Template file '#{@options.template_file}' is not valid: #{msg}"
    return false
  end #template_valid?()
  
  # Given an input get what its output should be based on the output opts
  # * +input+ - Input object, either file path, or maybe it is is @stdin
  # ==== Returns
  # * +@stdout+ - if @stdin
  # * All else build an output File using the prefix or the inplace extension
  def output_for_input(input)
    if @stdin == input
      return @stdout
    end
    if @options.prefix then
      return output_for_input_using_prefix(input)
    elsif @options.save_inplace then
      return output_for_input_using_inplace(input)
    end
    # if we got this far, we hopefully had only one input total, so go STDOUT
    if 1 == @options.input_files.size 
      return @stdout
    end
  end #output_for_input()
  
  # Get the path for the resulting output file using a prefix
  # * +input+ - the input pedigree file path
  # * Will use the @options.prefix to make the new file path
  def output_for_input_using_prefix(input)
    parts = base_dir_and_file_of(input)
    return File.join(parts.first,"#{@options.prefix}_#{parts.last}")
  end #output_for_input_using_prefix(input)
  
  # Get an output file using inplace replacement
  # * +input+ - input file path
  # * Will use the @options.backup_extension to make the new file path
  def output_for_input_using_inplace(input)
    if @options.backup_extension then
      return input + @options.backup_extension
    end
    return Tempfile.new('fix_vcf_generated_ped')
  end #output_for_input_using_inplace(input)
  
  # Get the base folder of the file
  def base_dir_and_file_of(input)
    [File.dirname(input),File.basename(input)]
  end #base_of(input)
  
  # Test the options given for output
  # * Either prefix, inplace, or neither
  # * Prefix requires input file(s), not STDIN
  # * Inplace require input file(s), not STDIN
  # * Inplace can have backup extension
  # * Neither prefix, nor inplace & only STDIN or only one file
  def output_options_valid?()
    if @options.prefix && @options.save_inplace then
      @stderr.puts "Selected output options are not valid, you can't have a prefix and save in place"
      return false      
    end
    
    if input_is_stdin?() && (@options.prefix || @options.save_inplace)
      @stderr.puts "Selected output options are not valid, you can't use STDIN with a prefix or saving inplace"
      return false      
    end
    
    if !(@options.prefix) && !(@options.save_inplace) && @options.input_files.size > 1 then
      @stderr.puts "Multiple inputs cannot go to STDOUT"
      return false
    end
    
    # by this point we either STDIN only, or Files with options that are valid
    return true
  end #output_options_valid
  
  # Test to see if the output path is valid
  def output_is_valid?(output)
    @stderr.puts "Testing #{output}" if @options.verbose 
    return true if @stdout == output
    if String == output.class then
      if File.exists?(output) then
        return (0 == File.size(output))
      else
        return true
      end
    elsif Tempfile == output.class
      return true
    end
    return false
  end #output_is_valid?(output)
  
  # Test if we are just doing STDING
  def input_is_stdin?()
    1 == @options.input_files.size && @stdin == @options.input_files[0]
  end #input_is_stdin?()
  
  # Just print the usage message
  def output_usage(out)
    out.puts <<-EOF
fix_vcf_generated_ped.rb -t PED_TEMPLATE -p OUTPUT_PREFIX INPUT_PEDIGREE(S)

Options:
 -h, --help             Display this help message
 -v, --version          Display the version information
 -V, --verbose          Increased verbosity of output
 -t, --template FILE    Specify the input template pedigree
 -a, --add_indivuals    Add any/all individuals in the template, regardless of them being present in original
 -p, --prefix STR       Specify the prefix the output file(s)
 -i, --inplace EXT      Save new file in place, saving original files with given extension.
                        Note original is not saved if no extension is given with -i option
                        
 Either use the prefix or inplace option for the new file, but not both
    
    EOF
  end

  def output_version(out)
    out.puts "#{File.basename(__FILE__)} Version: #{VERSION} Released: #{REVISION_DATE}"
  end
  
  def output_help(out)
    output_version(out)
    out.puts ""
    output_usage(out)
  end
  
end #PedFixerApp


if $0 == __FILE__
  PedFixerApp.new(ARGV.clone).run
end