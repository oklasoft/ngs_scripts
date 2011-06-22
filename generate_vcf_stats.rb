#!/usr/bin/env ruby1.9 -w -Ku
# encoding: UTF-8
#
# generate_vcf_stats.rb
# Created by Stuart Glenn on 2011-06-17T14:37:38-0500 
#
# == Synopsis
# Quick script to take soem number of VCFs and produce some stats for them.
# These stats include the vcf-stat output, a bed file, then some minor allele
# frequency information. We also intersect each with option UCSC style tracks
# of interest for the allele frequencying 
#
# == Inputs
#  - A file listing a name & pointer the the bed file for an interest region
#  - One or more VCFs
#
# == Usage
#  generate_vcf_stats.rb -t BED_TRACKS INPUT_VCF(S) 
#
#  For help use fix_vcf_generated_ped.rb -h
#
# == Options
#  -h, --help             Display this help message
#  -v, --version          Display the version information
#  -V, --verbose          Increased verbosity of output
#  -t, --tracks FILE      A file of track names & bed files to intersect the VCF
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

class VcfStatGeneratorApp
  VERSION       = "0.0.1-pre01"
  REVISION_DATE = "2011-06-17"
  AUTHOR        = "Stuart Glenn <Stuart-Glenn@omrf.org>"
  COPYRIGHT     = "Copyright (c) 2011 Oklahoma Medical Research Foundation"
  
  NUM_TRACK_FIELDS = 2
  TRACK_FILES_DELIMITER = /\s+/
  
  # Create an app object to run
  # * +args+ - Array of command line argument string
  # * +ios+ - Optional hash of :stdin, :stdout, :stderr to set to corresponding @std
  def initialize(args,ios = {})
    @args = args
    @error_message = ""
    set_inputs_outputs(ios)
    set_default_options()
  end

  # Run the application
  def run
    return false unless options_parsed?
    unless options_valid?
      output_usage(@stderr)
      exit(1)
    end

    @options.input_files.each do |vcf_path|
      unless process_vcf(vcf_path)
        @stderr.puts "Error processing #{vcf_path}: #{@error_message}"
        exit 1
      end
    end #each vcf
  end #run
  
  
  private
  
  # Does all the work for a vcf
  # Produces a basic bed file, makes intersections with tracks of interests,
  # calculated basaic minor allele frequency within & across those tracks,
  # produce a vcf-stats summary file as well
  # * *Args*    :
  #   - +vcf_path+ - The file system path to the VCF to process
  # * *Returns* :
  #   - false if it had problems, true if all good
  #   - It does also create 4 files in the same directory of the vcf_file
  # * *Raises* :
  #  - +Exception+ - Calls many sub methods which can screw up
  #
  def process_vcf(vcf_path)
    base_output_dir = output_base(vcf_path)
    return false if nil == base_output_dir
    
    new_file_prefix = vcf_file_name_prefix(vcf_path)
    return false if nil == new_file_prefix 

    return false unless make_bed(vcf_path,base_output_dir,new_file_prefix)
    
    return false unless produce_vcf_stats(vcf_path,base_output_dir,new_file_prefix)
    
    return false unless produce_allele_frequency_info(vcf_path,base_output_dir,new_file_prefix)

    return true
  end #process_vcf
  
  # Read the VCF and make a basic BED file in base_dir
  #
  # * *Args*    :
  #   - +vcf_path+ - Input VCF file path
  #   - +base_dir+ - The base folder into which to save our output
  #   - +new_prefix+ - The prefix to attach to our output bed file
  # * *Returns* :
  #   - true on success, false on errors
  # * *Raises* :
  #  - +Exception+ - Things might cause problems down in here
  #
  def make_bed(vcf_path,base_dir,new_prefix)
    @error_message = "Some trouble making the bed"
    File.open( File.join(base_dir,"#{new_prefix}.bed"), "w") do |out_bed|
      vcf_data_lines(vcf_path) do |data|
        vcf_data_to_bed_if_bed(data,out_bed)
      end
    end #bed file
    return true
  end #make_bed
  
  # Get lines of data from vcf into a hash & yield to the block
  #
  # * *Args*    :
  #   - +vcf_file+ - The VCF file
  #   - +block+ - The block which will be yieled with a data hash of the fields
  # * *Returns* :
  #   - Nothing of worth
  # * *Raises* :
  #  - +Exception+ -
  #
  def vcf_data_lines(vcf_file,&block)
    File.open(vcf_file).each do |line|
      next if line =~ /($^)|^#/
      data = {}
      (data[:chr],data[:pos],data[:id],data[:ref],data[:alt],data[:qual],
      data[:filter],data[:info]) = line.chomp.split(/\t/)
      block.call(data)
    end
  end #vcf_data_lines
  
  
  def vcf_data_to_bed_if_bed(vcf_data,out)
    return unless vcf_data && vcf_data[:chr] && vcf_data[:pos]
    result = []
    result << vcf_data[:chr].sub(/chr/i,'')
    result << vcf_data[:pos].to_i - 1
    result << "TODO"
    result << "#{vcf_data[:ref]}/#{vcf_data[:alt]}"
    #chr start end ref/alt
    out.puts result.join("\t")
  end
  
  # Call vcf-stats on the VCF saving those stats
  #
  # * *Args*    :
  #   - +vcf_path+ - Input VCF file path
  #   - +base_dir+ - The base folder into which to save our output
  #   - +new_prefix+ - The prefix to attach to our output vcf-stats file
  # * *Returns* :
  #   - true on success, false on errors
  # * *Raises* :
  #  - +Exception+ - Things might cause problems down in here
  #
  def produce_vcf_stats(vcf_path,base_dir,new_prefix)
    @error_message = "Some trouble calling vcf-stats"
    return false
  end #make_bed
  
  # Find some allele freqs & report them, as a whole & within tracks of interest
  #
  # * *Args*    :
  #   - +vcf_path+ - Input VCF file path
  #   - +base_dir+ - The base folder into which to save our output
  #   - +new_prefix+ - The prefix to attach to our output allele freq info files
  # * *Returns* :
  #   - true on success, false on errors
  # * *Raises* :
  #  - +Exception+ - Things might cause problems down in here
  #
  def produce_allele_frequency_info(vcf_path,base_dir,new_prefix)
    @error_message = "Some trouble calculating allele frequencies"
    return false
  end #make_bed
  
  
  # Given the path to a VCF file get what the prefix will be for all our new files
  #
  # * *Args*    :
  #   - +vcf_path+ - The path tring to the file being used
  # * *Returns* :
  #   - a string to use as the prefix on all newly generated files pertaining to this VCF, 
  #     or of there are troubles nil
  # * *Raises* :
  #  - ++ -
  #
  def vcf_file_name_prefix(vcf_path)
    @error_message = "Some sort of of trouble getting the prefix"
    file = File.basename(vcf_path,".vcf")
    if file =~ /\.recode$/
      file = File.basename(file,".recode")
    end
    return file
  end #vcf_file_name_prefix
  
  
  # Get the base directoy for output
  #
  # * *Args*    :
  #   - +vcf_path+ - The input VCF file
  # * *Returns* :
  #   - a directory path to save our output file for the given input vcf, nil if the path isn't good
  #
  def output_base(vcf_path)
    @error_message = "Unable to find/set base output directory"
    dir = File.dirname(vcf_path)
    unless File.writable?(dir)
      @error_message = "Base output dir, #{dir}, is not writeable"
      return nil
    end
    return dir
  end #output_base
  
  
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
      :tracks_file => nil,
      :verbose => false,
      :input_files => nil
    )
  end #set_default_options
  
  # Do the work of parsing out the options in the @arg array into the right
  # places into the @options
  # ==== Retuns
  # * +true+ - For good option parsing
  # * +false+ - On errors during parsing
  def options_parsed?
    opt_parser = OptionParser.new() do |opts|
      opts.on('-v','--version') { output_version(@stdout); exit(0) }
      opts.on('-h','--help') { output_help(@stdout); exit(0) }
      opts.on('-V', '--verbose')    { @options.verbose = true }

      opts.on("-t","--tracks", "=REQUIRED") do |tracks_file|
        @options.tracks_file = tracks_file
      end
    end

    opt_parser.parse!(@args) #rescue return false
    @options.input_files = @args unless @args.empty?
    return true
  end #options_parsed?

  # Make sure the options we have are correct & valid
  def options_valid?
    tracks_valid?() &&
    vcf_given?()
  end #options_valid?
  
  # Test to see that we were given at least one VCF argument
  #
  # * *Returns* :
  #   - true if we have at least one arg, false if not
  #
  def vcf_given?()
    unless @options.input_files && !@options.empty?
      @stderr.puts "Missing VCF(s) to process"
      return false
    end
    return true
  end #vcf_given?
  
  
  # Test the input tracks by checking it is readable & contains the valid 2 
  # columns, if given at all
  def tracks_valid?()
    return true unless @options.tracks_file
    msg = 'Unable to read file'
    if File.readable?(@options.tracks_file)
      msg = "Wrong format or no data"
      IO.foreach(@options.tracks_file) do |line|
        next if line =~ /^#/
        if NUM_TRACK_FIELDS == line.split(TRACK_FILES_DELIMITER).size
          return true
        else
          msg = "Wrong number of fields"
        end
      end
    end
    @stderr.puts "Interest Tracks file '#{@options.tracks_file}' is not valid: #{msg}"
    return false
  end #template_valid?()
  
  # Just print the usage message
  def output_usage(out)
    out.puts <<-EOF
#{File.basename($0)} -t BED_TRACKS INPUT_VCF(S)

Options:
 -h, --help             Display this help message
 -v, --version          Display the version information
 -V, --verbose          Increased verbosity of output
 -t, --tracks BEDS      Specify the input file with tracks of interest
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


  
end

if $0 == __FILE__
  VcfStatGeneratorApp.new(ARGV.clone).run
end