#!/usr/bin/env ruby1.9
#
# scatter_gather_gatk_variant_calling.rb
# Created by Stuart Glenn on 2011-06-07
#
# == Synopsis
# Quick script to wrap the work of running GATK scatter/gathter type || UnifiedGenotyper
#
# ==Author
# Stuart Glenn <Stuart-Glenn@omrf.org>
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

require 'ostruct'
require 'optparse'

class Interval
  attr_reader :chr, :start, :stop
  def initialize(chr,start,stop)
    @chr = chr
    @start = start.to_i
    @stop = stop.to_i
  end

  def size
    @stop - @start
  end

  def split
    mid = (((@stop-@start).abs)/2) + @start
    if mid <= @start || mid >= @stop || mid+1 >= @stop
      return self
    end
    return [self.class.new(@chr,@start,mid), self.class.new(@chr,mid+1,@stop)]
  end

  def to_s
    "#{@chr}:#{@start}-#{@stop}"
  end

  def self.intervals_from_file(file)
    res = []
    IO.foreach(file) do |l|
      res << self.new(*l.chomp.scan(/(.*):(\d+)-(\d+)/).first)
    end
    return res
  end
end

def normalize_interval_sizes(intervals,max_size)
  new_intervals = intervals.clone
  begin
    new_intervals = new_intervals.map do |i|
      if i.size > max_size then
        i.split
      else
        i
      end
    end.flatten
  end while new_intervals.detect {|i| i.size > max_size}
  new_intervals
end

def sliced_intervals()
  intervals = normalize_interval_sizes(Interval.intervals_from_file(@options.interval_list),@options.normalization_size)

  intervals_per_scatter = (intervals.size/@options.scatter_limit.to_f).ceil.to_i
  intervals.each_slice(intervals_per_scatter)
end

def output_version(out)
  out.puts "#{File.basename($0)} 1.0"
end

def output_help(out)
  out.puts "#{File.basename($0)} -o OUT_DIR -b BAM_LIST -i INTERVAL_LIST -r REF"
  out.puts ""
  out.puts <<-EOF
  -v, --version             Print version info
  -h, --help                Print this help
  -V, --verbose             Enable verbosity
  -D, --debug               Enable debuging
  -o, --output DIR          Save output in DIR
  -i, --interval LIST_FILE  File of intervals
  -b, --bam_list BAM_LIST   File listing the BAMs
  -r, --reference FASTQ     Fasta/q reference file against the BAMs were aligned
  -d, --dbsnp DBSNP         A VCF of SNP info for the UnifiedGenotyper
  -s, --scatter NUM         Scattering into NUM sub instances, defaults to 100
  --caller TYPE             Which type of caller to use, default UnifiedGenotyper, can be HaplotypeCaller
  --stand_call_conf VAL     Set stand_call_conf for GATK, defaults to 30.0
  --stand_emit_conf VAL     Set stand_emit_conf for GATK, defaults to 10.0
  -t, --threads NUM         Run the GATK with NUM threads, defaults to 1
  -n, --normalized NUM      Try to normalize the size intervals by be NUM size or less, defaults to 5000
EOF
end

def parse_opts(args)
  @options = OpenStruct.new(
    :output_base => Dir.pwd,
    :scatter_limit => 100,
    :stand_call_conf => 30.0,
    :stand_emit_conf => 10.0,
    :caller => 'UnifiedGenotyper',
    :verbose => false,
    :debug => false,
    :bam_list => nil,
    :interval_list => nil,
    :reference_path => nil,
    :snp_path => nil,
    :normalization_size => 5000,
    :threads => 1
  )
  
  opts = OptionParser.new() do |o|
    o.on('-v','--version') { output_version($stdout); exit(0) }
    o.on('-h','--help') { output_help($stdout); exit(0) }
    o.on('-V', '--verbose')    { @options.verbose = true }
    o.on('-D', '--debug')    { @options.debug = true }

    o.on("-o","--output", "=REQUIRED") do |output|
      @options.output_base = File.expand_path(output)
    end

    o.on("-s","--scatter", "=REQUIRED") do |scatter|
      @options.scatter_limit = scatter.to_i
    end

    o.on("-n","--normalized", "=REQUIRED") do |normalization_size|
      @options.normalization_size = normalization_size.to_i
    end

    o.on("-t","--threads", "=REQUIRED") do |threads|
      @options.threads = threads.to_i
    end

    o.on("--caller" "=OPTIONAL") do |conf|
      @options.caller = conf.to_s
    end

    o.on("--stand_call_conf", "=REQUIRED") do |conf|
      @options.stand_call_conf = conf.to_f
    end

    o.on("--stand_emit_conf", "=REQUIRED") do |conf|
      @options.stand_emit_conf = conf.to_f
    end

    o.on("-b","--bam_list", "=REQUIRED") do |bam|
      @options.bam_list = File.expand_path(bam)
    end

    o.on("-i","--interval", "=REQUIRED") do |interval|
      @options.interval_list = File.expand_path(interval)
    end

    o.on("-r","--reference", "=REQUIRED") do |path|
      @options.reference_path = File.expand_path(path)
    end

    o.on("-d","--dbsnp", "=REQUIRED") do |path|
      @options.snp_path = File.expand_path(path)
    end
  end

  opts.parse!(args)  

  unless @options.bam_list && @options.interval_list && @options.reference_path
    $stderr.puts "Need to at least give us a file listing bams & the master intervals & a refernce"
    $stderr.puts ""
    output_help($stderr)
    exit 1
  end  
end

parse_opts(ARGV)

STDERR.puts "#{@options}" if @options.verbose
name_base = File.split(@options.bam_list)[-1]

Dir.chdir(@options.output_base) do
  puts "Scatter/Gather Variants Starting #{$$}"
  work_dir = "#{Time.now.strftime("%Y%m%d")}_scatter_gather_work.#{$$}"
  raise "Can't make work dir" unless Dir.mkdir(work_dir)
  to_joins = []
  Dir.chdir(work_dir) do
    raise "Can't make log dir" unless Dir.mkdir("logs")
    
    threaded_queue = if @options.threads > 1
      "-pe threaded #{@options.threads}"
    else
      ""
    end
    
    sliced_intervals().each_with_index do |sliced_intervals,slice|
      slice = slice.to_s.rjust(@options.scatter_limit.to_s.length,"0")

      STDERR.puts "#{slice}\t#{sliced_intervals.join("; ")}" if @options.verbose

      sliced_interval_file = "sliced_interval_#{slice}.interval_list"
      File.open(sliced_interval_file,"w") do |f|
        f.puts sliced_intervals.join("\n")
      end

      output_file = File.join(Dir.pwd,"#{name_base}_#{slice}_variants.vcf")
      input_file = File.join(Dir.pwd,sliced_interval_file)
      
      snp_opt = case @options.snp_path 
        when nil
          ""
        when /\.rod$/
          "-D #{@options.snp_path}"
        when /\.vcf(\.gz)?$/
          "-D #{@options.snp_path} "
        else
          ""
      end

      cmd = "qsub #{threaded_queue} -p -10 -m e -o logs -b y -V -j y -cwd \
      -N #{name_base}_variants_#{slice} -l mem_free=4G,virtual_free=4G,h_vmem=6G \
gatk -T #{@options.caller} -glm BOTH -nt #{@options.threads} \
-A AlleleBalance \
-R #{@options.reference_path} #{snp_opt} \
-I #{@options.bam_list} \
-o #{output_file} \
-stand_call_conf #{@options.stand_call_conf} \
-stand_emit_conf #{@options.stand_emit_conf} \
-L #{input_file}"
      puts cmd
      system cmd
      sleep(1)

      to_joins << "-V #{output_file}"
    end #each slice
  end #work_dir
  
  slices_of_to_joins = to_joins.each_slice(100000/to_joins.first.length).to_a
  if 1 == slices_of_to_joins.size
    cmd = "qsub -m e -b y -V -j y -cwd -q all.q -N #{name_base}_variants_merge \\
  -l mem_free=4G,virtual_free=4G,h_vmem=6G gatk -T CombineVariants \\
  -genotypeMergeOptions UNSORTED \\
  -R #{@options.reference_path} \\
  -o #{name_base}_variants.vcf"
    cmd += " \\
#{to_joins.join(" \\\n")}
  "
    File.open("#{work_dir}-joiner.sh","w") do |f|
      f.puts "#!/usr/bin/env bash"
      f.puts "source /etc/profile.d/*.sh"
      f.puts "module load sge"
      f.puts "module load gatk/2.0-39-gd091f72"
      f.puts cmd
    end
  else
    alphabet = ("aa".."zz").to_a
    intermediate_vcfs_to_merge = []
    slices_of_to_joins.each_with_index do |slice_of_to_joins,index|
      cmd = "qsub -m e -b y -V -j y -cwd -q all.q -N #{alphabet[index]}_#{name_base}_variants_merge \\
      -l mem_free=4G,virtual_free=4G,h_vmem=6G gatk -T CombineVariants \\
      -genotypeMergeOptions UNSORTED \\
      -R #{@options.reference_path} \\
      -o #{alphabet[index]}_#{name_base}_variants.vcf"
      cmd += " \\
#{slice_of_to_joins.join(" \\\n")}
      "
      intermediate_vcf = File.join(Dir.pwd,"#{alphabet[index]}_#{name_base}_variants.vcf")
      intermediate_vcfs_to_merge << "-V #{intermediate_vcf}"
      File.open("#{work_dir}-#{alphabet[index]}-joiner.sh","w") do |f|
        f.puts "#!/usr/bin/env bash"
        f.puts "source /etc/profile.d/*.sh"
        f.puts "module load sge"
        f.puts "module load gatk/2.0-39-gd091f72"
        f.puts cmd
      end      
    end
    
    cmd = "qsub -m e -b y -V -j y -cwd -q all.q -N final_#{name_base}_variants_merge \\
    -l mem_free=4G,virtual_free=4G,h_vmem=6G gatk -T CombineVariants \\
    -genotypeMergeOptions UNSORTED \\
    -R #{@options.reference_path} \\
    -o #{name_base}_variants.vcf"
    cmd += " \\
#{intermediate_vcfs_to_merge.join(" \\\n")}
    "
    File.open("#{work_dir}-final-joiner.sh","w") do |f|
      f.puts "#!/usr/bin/env bash"
      f.puts "source /etc/profile.d/*.sh"
      f.puts "module load sge"
      f.puts "module load gatk/2.0-39-gd091f72"
      f.puts cmd
    end      
    
    
  end
end #base_dir
