#!/usr/bin/env ruby
#
#  Copyright (c) 2015 Stuart Glenn, Oklahoma Medical Research Foundation. (OMRF)
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
require 'tmpdir'
require 'json'

def set_default_opts
  @options = OpenStruct.new(
    :tmp_dir_base => Dir.pwd,
    :output_base => Dir.pwd,
    :verbose => false,
    :debug => false,
    :reference_index => nil,
    :reference_path => nil,
    :gtf_path => nil,
  )
end

def parse_args(args)
  opts = OptionParser.new() do |o|
    o.on('-v','--version') { output_version($stdout); exit(0) }
    o.on('-h','--help') { output_help($stdout); exit(0) }
    o.on('-V', '--verbose')    { @options.verbose = true }
    o.on('-D', '--debug')    { @options.debug = true }

    o.on("-o","--output", "=REQUIRED") do |output_destination|
      @options.output_base = File.expand_path(output_destination)
    end

    o.on("-t","--tmp", "=REQUIRED") do |topts|
      @options.tmp_dir_base = File.expand_path(topts)
    end

    o.on("-i","--index", "=REQUIRED") do |index_base|
      @options.reference_index = File.expand_path(index_base)
    end

    o.on("-r","--refrence", "=REQUIRED") do |ref|
      @options.reference_path = File.expand_path(ref)
    end

    o.on("-g","--gtf", "=REQUIRED") do |g|
      @options.gtf_path = File.expand_path(g)
    end

    o.on("-q","--qsub", "=REQUIRED") do |qopts|
      @options.qsub_opts = qopts
   end

  end

  opts.parse!(args) rescue return false
  @options.samples = @args
  return true
end

def main()
  set_default_opts()

  unless parse_args(ARGV)
    output_help($stderr)
    exit(1)
  end

  if ARGV.size <= 0
    output_help($stderr)
    exit(1)
  end

  index = (ENV['SGE_TASK_ID']||1).to_i - 1
  threads = (ENV['NSLOTS']) || 1

  lib = JSON.parse(ARGV.pop, {:symbolize_names=>true})[index]

  output = File.join(@options.output_base,
                     "#{index}.bam")

  return_val = 0
  tmp_prefix_dir = Dir.mktmpdir(index.to_s,@options.tmp_dir_base)
  begin
    Dir.chdir(tmp_prefix_dir) do
      puts "Aligning #{lib[:paths].join(",")}" if @options.verbose
      STDOUT.flush

      Dir.mkdir("1")
      Dir.chdir("1") do
        cmd = %W/STAR  --genomeDir #{@options.reference_index} --readFilesIn/
        cmd += lib[:paths]
        case File.extname(lib[:paths].first)
        when /\.gz\z/i
          cmd += ["--readFilesCommand", "zcat"]
        end
        cmd += %W/--runThreadN #{threads}/

        puts cmd.join(" ") if @options.verbose
        unless system(*cmd)
          $stderr.puts "Failure with STAR first pass"
          return_val = $?.exitstatus
        end
      end

      if 0 == return_val
        puts ""
        Dir.mkdir("index_pass2")
        Dir.chdir("index_pass2") do
          puts "Making new index using splice junction info from first pass" if @options.verrbose
          cmd = %W/STAR --runMode genomeGenerate --genomeDir . --genomeFastaFiles #{@options.reference_path}
                   --sjdbOverhang 99 --sjdbFileChrStartEnd ..\/1\/SJ.out.tab/
          if @options.gtf_path
            cmd += %W/--sjdbGTFfile #{@options.gtf_path}/
          end
          cmd += %W/--runThreadN #{threads}/
          puts cmd.join(" ") if @options.verbose
          unless system(*cmd)
            $stderr.puts "Failure with STAR genomeGenerate"
            return_val = $?.exitstatus
          end
          @options.star_index = Dir.pwd
        end
      end

      if 0 == return_val
        puts ""
        cmd = %W/STAR --outFileNamePrefix #{File.join(@options.output_base,index.to_s)}
                --outSAMtype BAM SortedByCoordinate --genomeDir #{@options.reference_index} --readFilesIn/
        cmd += lib[:paths]
        case File.extname(lib[:paths].first)
        when /\.gz\z/i
          cmd += ["--readFilesCommand", "zcat"]
        end
        cmd += %W/--runThreadN #{threads}/

        cmd += %W/--outStd BAM_SortedByCoordinate/
        picard = "picard AddOrReplaceReadGroups VALIDATION_STRINGENCY=LENIENT \
        MAX_RECORDS_IN_RAM=6000000 CREATE_INDEX=True COMPRESSION_LEVEL=8 \
        I=/dev/stdin O=#{output} SO=coordinate RGID=#{lib[:rg][:id]} RGLB=#{lib[:rg][:lb]} \
        RGPL=#{lib[:rg][:pl]} RGPU=#{lib[:rg][:pu]} RGSM=#{lib[:rg][:sm]}"

        cmd = cmd.join(" ") + "|" + picard

        cmd = "/bin/bash", "--norc", "--noprofile", "-v", "-o", "pipefail", "-o", "errexit", "-c", cmd
        puts cmd.join(" ") if @options.verbose
        unless system(*cmd)
          $stderr.puts "Failure with STAR second pass"
          return_val = $?.exitstatus
        end
      end
  end
  ensure
    FileUtils.remove_entry tmp_prefix_dir
    %w/_STARtmp Log.final.out Log.out Log.progress.out Log.std.out/.each do |f|
      file = "#{File.join(@options.output_base,index.to_s)}#{f}"
      FileUtils.remove_entry(file) if File.exists?(file)
    end
  end

  exit return_val
end

main()
