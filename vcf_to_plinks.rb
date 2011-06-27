#!/usr/bin/env ruby1.9 -w -Ku
# encoding: UTF-8
#
# vcf_to_plink.rb
# Created by Stuart Glenn on 2011-06-27T11:31:57-0500 
#
# == Synopsis
# Another hack of a script to take a VCF and make plink style ped/fam & also
# optionally make imputed versions of those based on additional data
#
# == Inputs
#  - A file listing the correct pedigree template structure
#  - A file of impute QC data
#  - One or more VCFs
#
# == Usage
#  vcf_to_plink.rb -t PEDIGREE_TEMPLATE -q QC_FILE -i INPUT_VCF(S) 
#
#  For help use fix_vcf_generated_ped.rb -h
#
# == Options
#  -h, --help             Display this help message
#  -v, --version          Display the version information
#  -V, --verbose          Increased verbosity of output
#  -t, --template FILE    A pedigree template file
#  -q, --quality FILE     A file of imputed QC data
#  -i, --impute           Also make an imputed ped/fam output
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

class VcfToPlink
  VERSION       = "0.0.1-pre01"
  REVISION_DATE = "2011-06-27"
  AUTHOR        = "Stuart Glenn <Stuart-Glenn@omrf.org>"
  COPYRIGHT     = "Copyright (c) 2011 Oklahoma Medical Research Foundation"

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

    return false unless create_plink_files_from_vcf(vcf_path,base_output_dir,new_file_prefix)
    
    if @options.do_impute then
      return false unless impute_vcf(vcf_path,output_base,new_file_prefix)
    end

    return true
  end #process_vcf
  
  # Creates (&fixes) a set of plink compatible ped & fam files from the VCF
  #
  # * *Args*    :
  #   - +vcf_path+ - The path to the input VCF file
  #   - +output_base+ - The base dir into which we will save our output
  #   - +new_file_prefix+ - The string we will prefix onto our new output files
  # * *Returns* :
  #   - true when things go well
  #   - false when things go wrong
  # * *Raises* :
  #  - +Exception+ -
  #
  def create_plink_files_from_vcf(vcf_path,output_base,new_file_prefix)
    @error_message = "Trouble making initial plink files"
    return make_plink(vcf_path,File.join(output_base,new_file_prefix)) &&
    fix_plink_ped(File.join(output_base,new_file_prefix))
  end #create_plink_files_from_vcf
  
  def fix_plink_ped(ped_prefix_path)
    cmd="fix_vcf_generated_ped.rb -i -t #{@options.template_pedigree} #{ped_prefix_path}.ped"
    if 0 == run_external_command(cmd,"fix_vcf_generated_ped failed")
      return true
    end
    return false
  end
  
  def make_plink(vcf_path,out_path)
    cmd= "vcftools --plink --vcf #{vcf_path} --out #{out_path} 1>/dev/null"
    if 0 == run_external_command(cmd,"vcftools to plink")
      File.unlink("#{out_path}.log")
      return true
    end
    return false
  end
  
  def run_external_command(cmd,msg="Command failed")
    unless system(cmd)
      @error_message = "#{msg}: #{$?}"
      return $?.exitstatus
    end
    return 0
  end
  
  
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
    @stdin = ios[:stdin] || $stdin
    @stdout = ios[:stdout] || $stdout
    @stderr = ios[:stderr] || $stderr
  end
  
  # Just make sure the app is happy with some steady & known default options
  # build the @options hash
  def set_default_options()
    @options = OpenStruct.new(
      :template_pedigree => nil,
      :verbose => false,
      :qc_file => nil,
      :do_impute => nil
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

      opts.on("-t","--template", "=REQUIRED") do |file|
        @options.template_pedigree = file
      end
  
      opts.on("-q","--quality", "=REQUIRED") do |file|
        @options.qc_file = file
      end

      opts.on("-i","--impute") do
        @options.do_input = true
      end
    end

    opt_parser.parse!(@args) #rescue return false
    @options.input_files = @args unless @args.empty?
    return true
  end #options_parsed?

  # Make sure the options we have are correct & valid
  def options_valid?
    template_pedigree_valid?() &&
    vcf_given?()
  end #options_valid?
  
  # Makes sure if we were given a pedigree template it is the right format
  #
  # * *Returns* :
  #   - true if no template OR if template is the right format
  #   - false if there is a template AND it is the wrong format
  #   - false if the given template file has any other problems (like can't be read)
  # * *Raises* :
  #  - +Exception+ -
  #
  def template_pedigree_valid?()
    return true if nil == @options.template_pedigree
    unless File.readable?(@options.template_pedigree)
      @stderr.puts "Template file '#{@options.template_pedigree}' is not actually readable"
      return false
    end
    IO.foreach(@options.tracks_file) do |line|
      next if line =~ /^#/
      if 7 == line.split(/\t/).size
        return true
      else
        @stderr.puts "Template file '#{@options.template_pedigree}' had the wrong number of fields"
        return false
      end
    end
    @stderr.puts "There was some other unknown problem with the template file '#{@options.template_pedigree}'"
    return false
  end #template_pedigree_valid
  
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
  
  # Just print the usage message
  def output_usage(out)
    out.puts <<-EOF
#{File.basename($0)} -t PEDIGREE_TEMPLATE -q QC_FILE -i INPUT_VCF(S)

Options:
 -h, --help             Display this help message
 -v, --version          Display the version information
 -V, --verbose          Increased verbosity of output
 -V, --verbose          Increased verbosity of output
 -t, --template FILE    A pedigree template file
 -q, --quality FILE     A file of imputed QC data
 -i, --impute           Also make an imputed ped/fam output
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
  VcfToPlink.new(ARGV.clone).run
end