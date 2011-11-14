module OMRF
  class VaastFeature
    def initialize()
      @variant_loci = {}
    end
    attr_accessor :id, :gene_name, :genome_permutation_p, :genome_permutation_95_lower, :genome_permutation_95_upper,
      :score, :rank
    attr_reader :variant_loci
    
    def blank?
      nil == @id && nil == @gene_name
    end
    
    def add_variant_locus_info(chrom,coord,type,data)
      @variant_loci[chrom] ||= {}
      @variant_loci[chrom][coord] ||= {}
      @variant_loci[chrom][coord][type] = data.clone
    end
    
    def add_data_from_vaast_output_line(line)
      md = line.match(/^>(.*)\s+(.*)/)
      if md
        self.id = md[1]
        self.gene_name = md[2]
      else
        add_data_line(line)
      end
    end
    
    private
    
    def add_locus_line(type,values)
      data = nil
      coord_index = nil
      case type
      when :br
        coord_index = 0
        data = {:change => values[1]}
      when :tr, :tu
        coord_index = 1
        data = {:score => values[0].to_f, :change => values[2]}
      end
      (coord,chrom) = values[coord_index].split("@")
      add_variant_locus_info(chrom,coord.to_i,type,data)
    end
    
    def add_data_line(line)
      return unless line =~ /^\w+.+:/
      (key,values) = line.split(":")
      return unless key
      values.strip!
      case key.downcase
        when "genome_permutation_0.95_ci"
          (lower,upper) = values.split(",")
          @genome_permutation_95_lower = lower.to_f
          @genome_permutation_95_upper = upper.to_f
        when "genome_permutation_p"
          @genome_permutation_p = values.to_f
        when "score"
          @score = values.to_f
        when "rank"
          @rank = values.to_f
        else
          add_locus_line(key.downcase.to_sym,values.split(/\t/))
      end
    end
  end #VaastFeature
  
  class VaastDatabase
    def find(chr,coordinate)
      return nil
    end
  end
  
  class VaastDatabaseOuputFiles < VaastDatabase
    def initialize(files)
      @vaast_data = {}
      files.each do |f|
        load_data_from_file(f)
      end
    end
    
    def size
      @vaast_data.size
    end
    
    def find(chr,coordinate)
      chrom_data = @vaast_data[chr]
      return nil unless chrom_data
      chrom_data[coordinate]
    end
    
    private
    
    def load_data_from_file(file)
      current_feature = VaastFeature.new()
      IO.foreach(file) do |line|
        next if line =~ /^##/
        if line =~ /^>/
          add_feature(current_feature)
          current_feature = VaastFeature.new()
        end
        current_feature.add_data_from_vaast_output_line(line.chomp)
      end #each line
      add_feature(current_feature)
    end #load_data_from_file
        
    def add_feature(feature)
      return if feature.blank?
      feature.variant_loci.each do |chrom,coords|
        @vaast_data[chrom] ||= {}
        coords.each do |coordinate,types|
          @vaast_data[chrom][coordinate] = 
          {
            :types => types,
            :id => feature.id,
            :gene_name => feature.gene_name, 
            :genome_permutation_p => feature.genome_permutation_p,
            :genome_permutation_95_lower => feature.genome_permutation_95_lower,
            :genome_permutation_95_upper => feature.genome_permutation_95_upper,
            :score => feature.score, 
            :rank => feature.rank
          }
        end
      end
    end
    
    
  end #VaastDatabaseOuputFiles
end #OMRF