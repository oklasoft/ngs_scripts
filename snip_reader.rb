# =======================================================================================================================
# expected file format:
# "lgs_code","Contig","Reference position","Consensus position","Variation type","Length","Reference","Variants","Allele variations","Frequencies","Counts","Coverage","Variant #1","Frequency of #1","Count of #1","Variant #2","Frequency of #2","Count of #2","Overlapping annotations","Amino acid change",
# =======================================================================================================================

require 'rubygems'
require 'fastercsv'
require 'optparse'

# -------------------------------------------------------------------
# mapping names from @regions_file to corresponding contig names
CHROMOSOME_TO_CONTIG = {
                        'chr1' => 'NC_000001 contig',
                        'chr2' => 'NC_000002 contig',
                        'chr3' => 'NC_000003 contig',
                        'chr4' => 'NC_000004 contig',
                        'chr5' => 'NC_000005 contig',
                        'chr6' => 'NC_000006 contig',
                        'chr7' => 'NC_000007 contig',
                        'chr8' => 'NC_000008 contig',
                        'chr9' => 'NC_000009 contig',
                        'chr10' => 'NC_000010 contig',
                        'chr11' => 'NC_000011 contig',
                        'chr12' => 'NC_000012 contig',
                        'chr13' => 'NC_000013 contig',
                        'chr14' => 'NC_000014 contig',
                        'chr15' => 'NC_000015 contig',
                        'chr16' => 'NC_000016 contig',
                        'chr17' => 'NC_000017 contig',
                        'chr18' => 'NC_000018 contig',
                        'chr19' => 'NC_000019 contig',
                        'chr20' => 'NC_000020 contig',
                        'chr21' => 'NC_000021 contig',
                        'chr22' => 'NC_000022 contig',
                        'chrX' => 'NC_000023 contig',
                        'chrY' => 'NC_000024 contig',
                        'chrM' => 'NC_001807 contig',
                       }
def generate_gene_data_from_regions_file
  # --------------------------------------------------------------------
  # create hash containing high and low contig values for each gene name
  # i.e.
  # @gene_data = {'NC_000001 contig' => {:low => 123456, :high => 234567, :name => 'PTPN22'},
  #              'NC_000002 contig' => {:low => 345678, :high => 456789, :name => 'CHROMEY'}}
  @gene_data = {}
  FasterCSV.read(@regions_file).each_with_index do |line, i|
    next if i == 0 || line[9].nil?
    chromosome = line[9]
    low = line[10].gsub(',', '').to_i
    high = line[11].gsub(',', '').to_i
    name = line[12]
    raise "Mapping for #{chromosome} not found; make sure regions file line #{i + 1} is properly formatted." unless CHROMOSOME_TO_CONTIG.keys.include?(chromosome)
    @gene_data[CHROMOSOME_TO_CONTIG[chromosome]] ||= []
    @gene_data[CHROMOSOME_TO_CONTIG[chromosome]] << {:low => low, :high => high, :name => name}
  end
end

def main
  parse_command_line_input
  generate_gene_data_from_regions_file
  contig_and_position = create_contig_and_position()
  results = create_results(contig_and_position)
  write_final_files(results)
end

def parse_command_line_input
  options = {}
  optparse = OptionParser.new do |opts|
    # set banner for help
    opts.banner = "Usage: snip_reader.rb -r regionsfile -b bedfile lgsfile1 lgsfile2 ..."

    opts.on( '-r', '--regions FILE', 'Specify name of regions file to use' ) do |file|
        options[:regions_file] = file
    end
    opts.on( '-b', '--bed FILENAME', 'Specify name to use for bed file; defaults to bed_file.bed' ) do |file|
        file << '.bed' unless file =~ /(\.bed)$/
        options[:bed_file] = file
    end
    opts.on( '-d', '--dip', 'Files to be read are dip files' ) do
        options[:dip] = true
    end
  end

  optparse.parse!
  @regions_file = options[:regions_file]
  @bed_filename = options[:bed_file] || 'bed_file.bed'
  @dip = options[:dip]
  @initial_files = if ARGV.size == 1
    # regex
    regexp = Regexp.new(ARGV.first)
    Dir.entries(File.dirname(__FILE__)).select { |filename| filename =~ regexp }
  else
    ARGV
  end

  raise "Must pass at least 1 lgs file and 1 regions file." unless @initial_files && @regions_file
  raise "Regions file must be a csv." unless @regions_file =~ /(\.csv)$/
end

def create_contig_and_position
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

  @lgs_codes = @initial_files.map do |file_name|
    file_name.split('_').first
  end

  @initial_files.each do |file_name|
    puts file_name
    lgs_code = file_name.split('_').first
    FasterCSV.read(file_name).each_with_index do |line, i|
      # skip header
      next if i == 0 || line[0].nil?

      contig = line[0].gsub(/\s.*/, '')
      ref_pos = line[1]
      cons_pos = line[2]
      var_type = line[3]
      length = line[4]
      reference = line[5]
      variants = line[6]
      allele_vars = line[7]
      frequencies = line[8]
      counts = line[9]
      coverage = line[10]
      variant_1 = line[11]
      frequency_of_1 = line[12]
      count_of_1 = line[13]
      variant_2 = line[14]
      frequency_of_2 = line[15]
      count_of_2 = line[16]
      overlapping_annotations = line[17]
      amino_acid_change = line[18]

      raise "unexpected formatting of line #{i + 1} in file '#{file_name}'" if [lgs_code, contig, ref_pos, reference, allele_vars, coverage, frequency_of_1].any?(&:nil?)

      contig_and_position[contig] ||= {}
      contig_and_position[contig][ref_pos] ||= []
      contig_and_position[contig][ref_pos] << {:lgs_code => lgs_code,
                                               :length => length,
                                               :reference => reference,
                                               :allele_vars => allele_vars,
                                               :coverage => coverage,
                                               :frequency_of_1 => frequency_of_1,
                                               :amino_acid_change => amino_acid_change
                                              }
    end
  end
  contig_and_position
end

def create_results(contig_and_position)
  results = []
  # check data and insert into results if it meets the criteria
  contig_and_position.each do |contig, positions_hash|
    positions_hash.each do |position, subjects|
      if subjects_meet_criteria?(subjects)
        results << {:normal_results => generate_array_of_results_data(contig, position, subjects), :pats_results => generate_array_of_pat_results_data(contig, position, subjects)}
      end
    end
  end
  results
end

def write_final_files(results)
  # write bed file
  File.open(@bed_filename, 'w') do |bed_file|
    # ---------------------------------------------------------------------------------------
    # write output files
    normal_results_heading = ["contig", "position", "number of subjects", "reference", "allele variations", "amino acid change", "gene name", "variation subjects"]
    pats_results_heading = ["contig", "position", "number of subjects", "reference", "allele variations", @lgs_codes].flatten
    FasterCSV.open("final_results.csv", "w") do |normal_csv|
      FasterCSV.open("pats_final_results.csv", "w") do |pats_csv|
        normal_csv << normal_results_heading
        pats_csv << pats_results_heading
        results.each do |hash|
          # length is last in line; leave out of output.
          # length is only necessary for bed file.
          # even then only necessary if @dip is true.
          normal_csv << hash[:normal_results][0..-2]
          pats_csv << hash[:pats_results]
          write_from_line_to_bed_file(hash[:normal_results], bed_file)
        end
      end
    end
    
  end
  
end

# =========================================================================================================================
# utility methods

def generate_array_of_results_data(contig, position, subjects)
  # 1. reference
  reference = subjects.first[:reference]
  # 2. allele variations
  allele_variations = subjects.map { |s| s[:allele_vars] }.uniq.join(',')
  # 3. variations and subjects belonging to them
  # ---------------------------------------------------------------------
  # create list of allele variations and subjects for each
  vars_to_subjects = {}
  subjects.each do |s|
    if s[:allele_vars] != ''
      vars_to_subjects[s[:allele_vars]] ||= []
      vars_to_subjects[s[:allele_vars]] << s[:lgs_code]
    end
  end
  vars_and_subjects = []
  vars_to_subjects.each do |allele_var, array_of_lgs_codes|
    lgs_codes = array_of_lgs_codes.join(',')
    vars_and_subjects << "#{allele_var}: #{lgs_codes}"
  end
  vars_and_subjects = vars_and_subjects.join('; ')
  # ---------------------------------------------------------------------
  # 4. amino acide change
  amino_acid_change = subjects.first[:amino_acid_change]
  # 5. gene name
  if @gene_data["#{contig}"]
    found_gene = @gene_data["#{contig}"].detect { |hash| (position.gsub(',','').to_i >= hash[:low] && position.gsub(',','').to_i <= hash[:high]) }
    gene_name = found_gene ? found_gene[:name] : ''
  else
    gene_name = ''
  end
  # 6. length
  # NOTE: this is only important for @dip files
  length = subjects.first[:length]
  # return pertinent information
  ["#{contig}", "#{position}", "#{subjects.count}", "#{reference}", "#{allele_variations}", "#{amino_acid_change}", "#{gene_name}", "#{vars_and_subjects}", "#{length}"]
end

def generate_array_of_pat_results_data(contig, position, subjects)
  # 1. reference
  reference = subjects.first[:reference]
  # 2. allele variations
  allele_variations = subjects.map { |s| s[:allele_vars] }.uniq.join(',')
  # 3. variations and subjects belonging to them
  # ---------------------------------------------------------------------
  # create list of allele variations and subjects for each
  vars_by_subject = @lgs_codes.map do |lgs_code|
    subject = subjects.detect { |sub| sub[:lgs_code] == lgs_code }
    if subject
      var = subject[:allele_vars] =~ /.*\/.*/ ? subject[:allele_vars].gsub('/', '') : subject[:allele_vars] + subject[:allele_vars]
    else
      var = reference + reference
    end
    var.split('').join("\t")
  end
  # ---------------------------------------------------------------------
  
  # return pertinent information
  ["#{contig}", "#{position}", "#{subjects.count}", "#{reference}", "#{allele_variations}", vars_by_subject].flatten
end

def write_from_line_to_bed_file(line, bed_file)
  contig = line[0]
  position = line[1].gsub(',','').to_i
  length = line[8].gsub(',', '').to_i
  number_of_subjects = line[2]
  reference = line[3]
  allele_variations = "(" + line[4].gsub(',', '_') + ")"
  gene_name = line[6] != '' ? line[6].gsub(' ', '_') : nil
  extras = [number_of_subjects, reference, allele_variations, gene_name].join('_').gsub(/_$/, '')
  # position output changes if using dip files
  if @dip
    bed_file.write("#{CHROMOSOME_TO_CONTIG[contig]}\t#{position}\t#{position + length}\t#{extras}\n")
  else
    bed_file.write("#{CHROMOSOME_TO_CONTIG[contig]}\t#{position - 1}\t#{position}\t#{extras}\n")
  end
end

def subjects_meet_criteria?(subjects)
  (subjects.size > 1) || (subjects.size == 1 && subjects.first[:frequency_of_1].to_f > 65 && subjects.first[:coverage].to_f > 10) || (subjects.size == 1 && subjects.first[:frequency_of_1].to_f < 65 && subjects.first[:coverage].to_f > 20)
end
# =============================================================================================================
# EXECUTE PROGRAM

main()
