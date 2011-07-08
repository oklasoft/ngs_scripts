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
#  - A file listing loci info of the target VCF(s): chr start stop name
#  - A file listing associated genetic map file locations for each chromosome: chr path_to_file
#  - A file listing the correct template pedigree file for the VCF made ped file
#  - A base directory contiaing genotype data for impute, named by locus (LOCUS.gen)
#  - A base directory contiaing sample data for impute, named by locus (LOCUS.sample)
#    This base dir will also be used to find a pedigree template (LOCUS.ped) to fix the imputed pedigree
#  - One or more VCFs
#
# == Usage
#  vcf_to_plink.rb -t PEDIGREE_TEMPLATE -q QC_FILE -i INPUT_VCF(S) 
#
#  For help use vcf_to_plink.rb -h
#
# == Options
#  -h, --help             Display this help message
#  -v, --version          Display the version information
#  -V, --verbose          Increased verbosity of output
#  -t, --template FILE    A pedigree template file
#  -g, --genotype DIR     A directory containing the genotypes named by loci for imputation
#  -i, --impute           Also make an imputed ped/fam output
#  -l, --loci FILE        A file with information about possible loci
#  -m, --map FILE         A file matching chromosome with associated path to a map file
#  -s, --sample DIR       A directory containing the sample info files named by loci for gtool
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
    
    unless load_loci()
      @stderr.puts "Error loading loci info: #{@error_message}"
      exit 1
    end

    unless load_genetic_maps()
      @stderr.puts "Error loading mapping info: #{@error_message}"
      exit 1
    end

    had_errors = false

    @options.input_files.each do |vcf_path|
      unless process_vcf(vcf_path)
        @stderr.puts "Error processing #{vcf_path}: #{@error_message}"
        had_errors = true
      end
    end #each vcf
    return !had_errors
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
      return false unless impute_vcf(vcf_path,base_output_dir,new_file_prefix)
    end

    return true
  end #process_vcf
  
  # Reads in the loci file to a hash used later to find matches for
  # the VCF files based on possible matches between file names & loci names
  #
  # * *Returns* :
  #   - true on no problems
  #   - false on problems
  #   - Also creates & fills the @loci hash
  # * *Raises* :
  #  - +Exception+ -
  #
  def load_loci()
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
  
  # Read in the paths to the genetic maps save to hash keyed by chr
  #
  # * *Returns* :
  #   - true on success & @genetic_maps is filled
  #   - false
  #
  def load_genetic_maps()
    @genetic_maps = {}
    IO.foreach(@options.genetic_map_file) do |line|
      next if line =~ /^#/ || line =~ /^$/
      (chr,path) = line.chomp.split(/\s+/)
      @genetic_maps[chr] = path
    end
    return true
  rescue => error
    @error_message = error.message
    return false
  end #load_genetic_maps
  
  
  # Imputes the given VCF & fixes up some pential problems
  #
  # * *Args*    :
  #   - +vcf_path+ - The path to the input VCF file
  #   - +output_base+ - The base dir into which we will save our output
  #   - +new_file_prefix+ - The string we will prefix onto our new output files
  # * *Returns* :
  #   - true when things go well
  #   - false when they go porly
  # * *Raises* :
  #  - +Exception+ -
  #
  def impute_vcf(vcf_path,output_base,new_file_prefix)
    out_path = File.join(output_base,new_file_prefix)
    
    locus = find_locus(new_file_prefix)
    unless locus
      @error_message = "We were unable to find a possible locus for #{new_file_prefix} from the loci file"
      return false
    end
    
    hapmap_chromsome_panel_file = @genetic_maps[locus[:chr]]

    population_locus_genotype_file = File.join(@options.genotype_dir,"#{locus[:name]}.gen")
    unless check_file(population_locus_genotype_file,{:readable => true})
      @error_message = "The genotype file we tried to use (#{population_locus_genotype_file}) is not readable"
      return false
    end
    
    gtool_sample_file = File.join(@options.samples_dir,"#{locus[:name]}.sample")
    imputed_pedigree_template = File.join(@options.samples_dir,"#{locus[:name]}.ped")
    unless check_file(gtool_sample_file,{:readable => true})
      @error_message = "The sample file we tried to use (#{gtool_sample_file}) is not readable"
      return false
    end
    
    if create_impute_input(vcf_path,out_path) && convert_indels_to_ac_snp_in_impute_legend(out_path)
      if run_impute(out_path,locus[:start],locus[:stop],hapmap_chromsome_panel_file,population_locus_genotype_file)
        return convert_impute_ouput_to_plink(out_path,gtool_sample_file) &&
          fix_imputed_ped_file("#{out_path}_imputed.ped",imputed_pedigree_template) &&
          fix_imputed_map_file("#{out_path}_imputed.map",locus[:chr])
      else
        @stderr.puts "WARN: #{vcf_path} failed to impute"
        clean_impute_files(out_path)
        %w/_impute_legend_indel_changes.txt _imputed_warnings/.each do |ext|
          File.unlink("#{out_path}#{ext}") if File.exists?("#{out_path}#{ext}")
        end
        return true
      end
    end
    return false
  end #impute_vcf
  
  def clean_impute_files(input_prefix)
    %w/_imputed _imputed_info _imputed_info_by_sample .impute.hap .impute.hap.indv .impute.legend/.each do |ext|
      File.unlink("#{input_prefix}#{ext}") if File.exists?("#{input_prefix}#{ext}")
    end
  end
  
  # Attempt to match the filename against our @loci
  #
  # * *Args*    :
  #   - +filename+ - A string we use to find possible mathes against
  # * *Returns* :
  #   - the matching locus info if found
  #   - nil if no match
  #
  def find_locus(filename)
    @loci.each do |locus_name,locus_data|
      if filename =~ /^#{locus_name}$/
        return locus_data
      end
    end
    return nil
  end #find_locus
  
  
  # Fix the map by readding chromoroms & renaming "----" "SNPS"
  #
  # * *Args*    :
  #   - +map_file+ - The file we'll change
  #   - +chr+ - The current chromosome
  # * *Returns* :
  #   - True when good, false when not
  # * *Raises* :
  #  - +Exception+ - real helpful I know
  #
  def fix_imputed_map_file(map_file,chr)
    @error_message = "Trouble fixing imputed map file '#{map_file}'"
    tmp_file = "#{map_file}.#{$$}"
    unless check_file(tmp_file,{:exists => false,:writable => true})
      @error_message = "The temp map file we would use can't be created"
      return false
    end
    File.open(tmp_file,"w") do |output|
      IO.foreach(map_file) do |line|
        fix_imputed_map_line_to(line.chomp,chr,output)
      end
    end
    File.unlink(map_file)
    File.rename(tmp_file,map_file)
    return true
  end #fix_imputed_map_file

  # set the chr & put position for rs if not there
  def fix_imputed_map_line_to(line,new_chr,output)
    return if line =~ /^$/ || line =~ /^#/
    parts = line.split(/\t/)
    parts[0] = new_chr
    parts[1] = parts[3] if "---" == parts[1]
    output.puts parts.join("\t")
  end
  
  
  # Fix that newly imputed ped file. Fixing is defined as:
  # 1- Chancing the N N unknown genotypes to 0 0
  # 2- Correcting the pedigree with full/correct info from template
  #
  # * *Args*    :
  #   - +ped_file+ - The input pedigree file to fix
  #   - +pedigree_template+ - The template pedigree file used to correct families
  # * *Returns* :
  #   - true on success
  #   - false if a scumbag caused problems
  # * *Raises* :
  #  - +Exception+ - Yea, real generic
  #
  def fix_imputed_ped_file(ped_file,impute_pedigree_template = nil)
    @error_message = "Trouble fixing imputed ped file '#{ped_file}'"
    tmp_file = "#{ped_file}.#{$$}"
    unless check_file(tmp_file,{:exists => false,:writable => true})
      @error_message = "The temp ped file we would use can't be created"
      return false
    end
    File.open(tmp_file,"w") do |output|
      IO.foreach(ped_file) do |line|
        fix_imputed_pedigree_line_to(line.chomp,output,nil != impute_pedigree_template)
      end
    end
    File.unlink(ped_file)
    File.rename(tmp_file,ped_file)
    if impute_pedigree_template
      cmd="fix_vcf_generated_ped.rb -i .backup -t #{impute_pedigree_template} #{ped_file}"
      if 0 == run_external_command(cmd,"fix_vcf_generated_ped for imputed ped failed (#{cmd})")
        return true
      end
      return false
    end
    return true
  end #fix_imputed_ped_file


  def fix_imputed_pedigree_line_to(line,output,put_indiv_id_as_ped_id=false)
    (ped_id,indiv_id,mother_id,father_id,sex,phenotype,*alleles) = line.split(/\t/)
    alleles.map! { |x| "N N" == x ? "0 0" : x}
    if put_indiv_id_as_ped_id
      ped_id = indiv_id
    end
    output.puts "#{ped_id}\t#{indiv_id}\t#{mother_id}\t#{father_id}\t#{sex}\t#{phenotype}\t#{alleles.join("\t")}"
  end
  
  # Convert the output of impute using gtool back to plink style files
  #
  # * *Args*    :
  #   - +input_prefix+ - What our output prefix setting to impute was, used to find our files
  # * *Returns* :
  #   - true on success
  #   - false on failures
  # * *Raises* :
  #  - +Exception+ -
  #
  def convert_impute_ouput_to_plink(input_prefix,input_sample_file)
    @error_message = "Trouble converting impute output to plink file '#{__LINE__}'"
    new_prefix = "#{input_prefix}_imputed"
    cmd=<<-EOF
gtool --log /dev/null -G --g #{new_prefix} --ped #{new_prefix}.ped --map #{new_prefix}.map \
--phenotype Phenotype \
--threshold 0.8 \
--s #{input_sample_file} 2>/dev/null 1>/dev/null
EOF
    if 0 == run_external_command(cmd,"gtool")
      # clean up impute files
      clean_impute_files(input_prefix)
      return true
    end
    return false
  end #convert_impute_ouput_to_plink
  
  
  # Do the bit of calling out & running impute
  #
  # * *Args*    :
  #   - +input_prefix+ - The base prefix to the input files from VCF, should be
  #                      Able to append .impute.hap & .impute.legend
  #   - +locus_start+ - Position the locus starts
  #   - +locus_stop+ - End position of this locus
  #   - +hapmap_chromsome_panel_file+- Path to the hampmap file for this impute
  #   - +population_locus_genotype_file+- File for some genotypes
  # * *Returns* :
  #   - true when things go well
  #   - false when they go poorly
  # * *Raises* :
  #  - ++ -
  #
  def run_impute(input_prefix,locus_start,locus_stop,hapmap_chromsome_panel_file,
                 population_locus_genotype_file)
    @error_message = "Trouble running impute '#{__LINE__}'"
    cmd=<<-EOF
impute -h #{input_prefix}.impute.hap -l #{input_prefix}.impute.legend \
-fix_strand_g \
-o #{input_prefix}_imputed \
-int #{locus_start} #{locus_stop} \
-m #{hapmap_chromsome_panel_file} \
-g #{population_locus_genotype_file} 2>/dev/null 1>/dev/null
EOF
    if 0 == run_external_command(cmd,"impute")
      return true
    end
    return false
  end #run_impute
  
  
  
  # Perform the checks as specified by modes
  #
  # * *Args*    :
  #   - +path+ - The path to the file
  #   - +modes+ - A hash of modes that need checking
  # * *Returns* :
  #   - True if the modes pass, false otherwise
  # * *Raises* :
  #  - +Exception+ -
  #
  def check_file(path,modes={})
    passing = true
    modes.each do |check,pass_fail|
      case check
        when :exists
          passing = passing & pass_fail ? File.exists?(path) : !File.exists?(path)
        when :readable
          passing = passing &  pass_fail ? File.readable?(path) : !File.readable?(path)
        when :writable
          if !File.exists?(path)
            passing = passing & pass_fail ? File.writable?(File.dirname(path)) : !File.writable?(File.dirname(path))
          else
            passing = passing & pass_fail ? File.writable?(path) : !File.writable?(path)            
          end
        else
          raise "Unsupported check: #{check}"
      end
    end
    return passing
  end #check_file
  
  # We redo the .impute.legend file to remove any indels
  # They are replaced by an A->C SNP
  # We also create a log file denoting each change we made
  #
  # * *Args*    :
  #   - +input_prefix+ - The prefix of the file to read, we'll append .impute.legend
  #                      Also used for naming the new file .impute.legend.tmp & log
  # * *Returns* :
  #   - true when things go well
  #   - false when things go wrong
  # * *Raises* :
  #  - +Exception+ -
  #
  def convert_indels_to_ac_snp_in_impute_legend(input_prefix)
    @error_message = "Trouble fixing imputed indels '#{__LINE__}'"
    input_file = "#{input_prefix}.impute.legend"
    unless check_file(input_file,{:readable => true})
      @error_message = "The input legend file, #{input_file}, is not readable"
      return false
    end
    
    new_file = "#{input_prefix}.impute.legend.without_indels.#{$$}"
    unless check_file(new_file,{:exists => false,:writable => true})
      @error_message = "We will be unable to make our new legend file: #{new_file}"
      return false
    end
    
    log_file = "#{input_prefix}_impute_legend_indel_changes.txt"
    unless check_file(log_file,{:exists => false,:writable => true})
      @error_message = "We will be unable to log our changes to: #{log_file}"
      return false
    end
    
    File.open(new_file,"w") do |output|
      output.puts "ID pos allele0 allele1"
      File.open(log_file,"w") do |log|
        IO.foreach(input_file) do |line|
          translate_legend_line_without_indels_with_log(line.chomp,output,log)
        end
      end
    end
    File.unlink(input_file)
    File.rename(new_file,input_file)
    return true
  end #convert_indels_to_ac_snp_in_impute_legend

  def translate_legend_line_without_indels_with_log(line,output,log)
    return if line =~ /^$/
    return if line =~ /^ID/
    (id,pos,ref,variant) = line.split(/\s+/)
    if ref.length > 1 || variant.length > 1
      log.puts [id,pos,ref,variant].join("\t")
      ref = "A"
      variant = "C"
    end
    output.puts "#{id} #{pos} #{ref} #{variant}"
  end
  
  
  # Calls out to vcftools to make the initial impute input files from the VCF
  #
  # * *Args*    :
  #   - +vcf_path+ - The path to the input VCF file
  #   - +output_path+ - The new file output path prefix stuff
  # * *Returns* :
  #   - true when things go well
  #   - false when things go wrong
  # * *Raises* :
  #  - +Exception+ -
  #
  def create_impute_input(vcf_path,out_path)
    cmd="vcftools --IMPUTE --vcf #{vcf_path} --out #{out_path} 1>/dev/null"
    if 0 == run_external_command(cmd,"vcftools to impute")
      File.unlink("#{out_path}.log")
      return true
    end
    return false
  end #create_impute_input
  
  
  
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
    return true unless @options.template_pedigree
    cmd="fix_vcf_generated_ped.rb -i .backup -t #{@options.template_pedigree} #{ped_prefix_path}.ped"
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
      :genotype_dir => nil,
      :samples_dir => nil,
      :do_impute => nil,
      :loci_file => nil,
      :genetic_map_file => nil
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

      opts.on("-l","--loci", "=REQUIRED") do |file|
        @options.loci_file = file
      end

      opts.on("-m","--maps", "=REQUIRED") do |file|
        @options.genetic_map_file = file
      end
  
      opts.on("-g","--genotype", "=REQUIRED") do |file|
        @options.genotype_dir = file
      end

      opts.on("-s","--sample", "=REQUIRED") do |file|
        @options.samples_dir = file
      end

      opts.on("-i","--impute") do
        @options.do_impute = true
      end
    end

    opt_parser.parse!(@args) #rescue return false
    @options.input_files = @args unless @args.empty?
    return true
  end #options_parsed?

  # Make sure the options we have are correct & valid
  def options_valid?
    template_pedigree_valid?() &&
    required_data_for_impute_present &&
    vcf_given?()
  end #options_valid?
  
  # Check if we are imputing we have the right additional options
  #
  # * *Returns* :
  #   - true when all the opts are present & good
  #   - false if there are missing or problematic opts
  # * *Raises* :
  #  - +Exception+ -
  #
  def required_data_for_impute_present()
    return true unless @options.do_impute

    unless @options.loci_file && check_file(@options.loci_file,{:readable => true})
      @stderr.puts "We need a loci file to know the genomic interval to impute"
      return false
    end
    
    unless @options.genetic_map_file && check_file(@options.genetic_map_file,{:readable => true})
      @stderr.puts "We need a listing of genetic map files to have maps to impute"
      return false
    end
    
    unless @options.genotype_dir
      @stderr.puts "We need a directory of genotype files"
      return false
    end
    
    unless @options.samples_dir
      @stderr.puts "We need a directory of sample files"
      return false
    end
    
    return true
  end #required_data_for_impute_present
  
  
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
    IO.foreach(@options.template_pedigree) do |line|
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
#{File.basename($0)} -t PEDIGREE_TEMPLATE [-i] INPUT_VCF(S)

Options:
 -h, --help             Display this help message
 -v, --version          Display the version information
 -V, --verbose          Increased verbosity of output
 -V, --verbose          Increased verbosity of output
 -t, --template FILE    A pedigree template file for the VCF generated ped file
 -g, --genotype DIR     A directory containing the genotypes for imputation
 -i, --impute           Also make an imputed ped/fam output
 -l, --loci FILE        Read chromsome & interval info for the loci from FILE. Required for imputing
 -m, --map FILE         Read chromsome & associated map file from FILE. Required for imputing
 -s, --sample DIR       A directory containing the sample info files named by loci for gtool
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
  if VcfToPlink.new(ARGV.clone).run
    exit 0
  else
    $stderr.puts "Errors in running, check output"
    exit 1
  end
end