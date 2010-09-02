#!/usr/bin/env ruby

unless Array.public_methods.include?("none?")
  class Array
  
    def none?(&block)
      !any?(&block)
    end
  
  end
end

# =======================================================================================================================
# expected file format:
# "lgs_code","Contig","Reference position","Consensus position","Variation type","Length","Reference","Variants","Allele variations","Frequencies","Counts","Coverage","Variant #1","Frequency of #1","Count of #1","Variant #2","Frequency of #2","Count of #2","Overlapping annotations","Amino acid change",
# =======================================================================================================================

require 'rubygems'
require 'fastercsv'
require 'optparse'

# -----------------------------------------------------------------
# Parse Options
# hash to hold options parsed from commandline by OptionParser
options = {}
optparse = OptionParser.new do |opts|
  # set banner for help
  opts.banner = "Usage: snip_reader.rb -n lgsNOTfile1 -n lgsNOTfile2 lgsANDfile1 lgsANDfile2 ..."
  
  opts.on( '-n', '--not FILE', 'Specify name of file to leave out of intersection. Can be used multiple times.' ) do |file|
      options[:not_files] ||= []
      options[:not_files] << file
  end
  opts.on( '-o', '--out NAME', 'Specify name of output directory; defaults to the specified intersection.') do |name|
      options[:outdirectory_name] = name
  end
  
end

optparse.parse!
NOT_FILES = options[:not_files] || []
INTERSECT_FILES = ARGV

raise "Must pass 2 or more intersection files" unless INTERSECT_FILES.size >= 2
raise "All filenames must begin with the subject's lgscode followed by '_'" unless (INTERSECT_FILES + NOT_FILES).all? { |filename| filename =~ /^lgs\d+_/ }

COUNTING_CODES = INTERSECT_FILES.map { |filename| filename.split('_').first }
NON_COUNTING_CODES = NOT_FILES.map { |filename| filename.split('_').first }

# display intersection info.
INTERSECTION_NAME = COUNTING_CODES.join(" AND ") + NON_COUNTING_CODES.map { |not_file| " NOT #{not_file}" }.join
OUTDIRECTORY_NAME = options[:outdirectory_name] || INTERSECTION_NAME
puts INTERSECTION_NAME

# -------------------------------------------------------------------
# hash tracking each contig, its positions, and the subjects in each position
# i.e.
# contig_and_position = {'NC_000001' => {'123,456' => [{:lgs_code => 'lgs000600', :reference => 'A',...},
#                                                      {:lgs_code => 'lgs002324', :reference => 'G',...}
#                                                     ],
#                                        '234,567' => [{:lgs_code => 'lgs123456', :reference => 'C',...}]
#                                       },
#                        'NC_000002' => {'456,678' => [{:lgs_code => 'lgs356456', :reference => 'T',...}]}
#                       }
# -------------------------------------------------------------------

contig_and_position = {}

(INTERSECT_FILES + NOT_FILES).each do |file_name|
  # extract lgscode from file_name
  raise "file #{file_name} must being with lgscode_" unless file_name =~ /^lgs\d+_/
  lgscode = file_name.split('_').first
  
  FasterCSV.read(file_name).each_with_index do |line, i|
    # skip header
    next if i == 0 || line[1].nil?
    
    contig, ref_pos, cons_pos, var_type, length, reference, variants, allele_vars, frequencies, counts, coverage, variant_1, frequency_of_1, count_of_1, variant_2, frequency_of_2, count_of_2, overlapping_annotations, amino_acid_change = line
    
    raise "unexpected formatting of line #{i + 1} in file '#{file_name}'" if (contig.nil? || ref_pos.nil? || reference.nil? || allele_vars.nil? || coverage.nil? || frequency_of_1.nil? )
    
    contig = contig.gsub('chr0', 'chr').gsub(' contig', '')
    ref_pos = ref_pos.gsub(",", '')
    
    contig_and_position[contig] ||= {}
    contig_and_position[contig][ref_pos] ||= []
    contig_and_position[contig][ref_pos] << {:lgscode => lgscode,
                                             :length => length,
                                             :reference => reference,
                                             :allele_vars => allele_vars,
                                             :coverage => coverage,
                                             :frequency_of_1 => frequency_of_1,
                                             :overlapping_annotations => overlapping_annotations,
                                             :amino_acid_change => amino_acid_change
                                            }
  end
end
# should not change any after this
contig_and_position.freeze

# -------------------------------------------------------------------
# detect and save positions that are common to all files

results_1 = []
results_2 = []
results_3 = []

contig_and_position.each do |contig, hash_by_position|
  hash_by_position.each do |position, subjects|
    lgscodes = subjects.map { |subject| subject[:lgscode] }
    if(COUNTING_CODES.all? { |code| lgscodes.include?(code) } && NON_COUNTING_CODES.none? { |code| lgscodes.include?(code) })
      reference_allele = subjects.first[:reference]
      # variants sometimes look whacky, i.e. 'A/T'
      # if one of the values is equal to the reference_allele, exclude it.
      # if neither is equal to the reference_allele, include it as is.
      crazy_variants = subjects.map { |subject| subject[:allele_vars] }.uniq
      variants = crazy_variants.map do |var|
        var1, var2 = var.split("/")
        if var2
          if var1 != reference_allele && var2 != reference_allele
            var
          else
            [var1, var2].detect { |value| value != reference_allele }
          end
        else
          var1
        end
      end.flatten.uniq
      overlapping_annotations = subjects.map { |subject| subject[:overlapping_annotations] }.uniq.join("; ")
      amino_acid_change = subjects.map { |subject| subject[:amino_acid_change] }.uniq.join("; ")
      variants.each do |variant|
        results_1 << [contig, position, reference_allele, variant, overlapping_annotations, amino_acid_change]
        results_2 << [contig, position, '1', reference_allele + "/" + variant]
        results_3 << [contig, (position.to_i - 1).to_s, position]
      end
    end
  end
end

# ---------------------------------------------------------------------------------------
# write output file
file1_heading = ["contig", "position", "reference", "allele variations", "overlapping annotations", "amino acid change"]
file2_heading = ["contig", "position", "ref/var"]
file3_heading = ["contig", "position-1", "position"]

Dir.mkdir(OUTDIRECTORY_NAME)
Dir.open(OUTDIRECTORY_NAME) do |dir|
  FasterCSV.open(dir.path + "/for_graham.csv", "w") do |csv|
    csv << file1_heading
    results_1.each { |line| csv << line }
  end
  FasterCSV.open(dir.path + "/for_program.csv", "w") do |csv|
    csv << file2_heading
    results_2.each { |line| csv << line }
  end
  File.open(dir.path + "/bed_file.bed", "w") do |file|
    file.write(file3_heading.join("\t") + "\n")
    results_3.each { |line| file.write(line.join("\t") + "\n") }
  end
  File.open(dir.path + "/intersection_info.txt", "w") do |file|
    file.write(INTERSECTION_NAME + "\n")
  end
end
