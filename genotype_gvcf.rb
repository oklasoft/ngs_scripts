#!/usr/bin/env ruby
# Copyright (c) 2015 Stuart Glenn, Oklahoma Medical Research Foundation
# Licensed under 3 Clause BSD
# Wrapper around various GATK commands to get VCF output from some g.vcf
# Basic workflow is
# 1- GenotypeGVCFs one or more g.vcfs (listing in .list or separately)
# 2- Split SNPs & Indels apart & VQSR each separately
# 3- Remerge those VQSRed sets
# 4- Recode with a minDP

require 'yaml'
require 'optparse'
require 'thread'
require 'fileutils'

STEPS_DIRS_FILES = {
  genotypegvcf:{dir:"00_original_vcf", files:[], ext:lambda {|f|["#{f}-raw.vcf.gz"]}},
  split:{
    dir:"01_split_snp_indel",
    ext:lambda { |f| ["#{f}-snps.vcf.gz","#{f}-indels.vcf.gz"] },
    files:[],
    types:%w/SNP INDEL/},
  vqsr:{
    dir:"01_split_snp_indel",
    ext:lambda { |f| "#{f}-vqsr_recalibrated.vcf" },
    files:[],
    types:%w/snp indel/},
  merge:{dir:"02_remerged_vcf", files:[], ext:lambda {|f|["#{f}-vqsr-recalibrated.vcf.gz"]}},
  recode:{dir:"03_recoded", files:[], ext:lambda {|f,d|["#{f}-removed-mindp#{d}-recoded.vcf"]}},
  tabix:{dir:"03_recoded", files:[], ext:lambda {|f|["#{f}.vcf.gz"]}},
}

def which(cmd)
  exts = ENV['PATHEXT'] ? ENV['PATHEXT'].split(';') : ['']
  ENV['PATH'].split(File::PATH_SEPARATOR).each do |path|
    exts.each { |ext|
      exe = File.join(path, "#{cmd}#{ext}")
      return exe if File.executable?(exe) && !File.directory?(exe)
    }
  end
  return nil
end

def parsed_opts(method=:parse)
  options = {
    debug:false,
    verbose:false,
    dp:10,
    mem:90,
    threads:6,
    output_base:"variants",
    inputs:[]
  }

  op = OptionParser.new do |o|
    o.banner = "Usage: #{File.basename(__FILE__)} -a CONFIG_YAML -i INPUT"
    o.on("-c","--config YAML","Specify project config") do |c|
      options[:conf_file] = File.expand_path c
    end

    o.on("-m",'--memory GB',OptionParser::DecimalInteger,"Specify GB memory limit for GATK, default #{options[:mem]}") do |m|
      options[:mem] = m
    end

    o.on("-t","--threads NUM",OptionParser::DecimalInteger,"Run GATK GenotypeGVCFs with NUM threads, default #{options[:threads]}") do |t|
      options[:threads] = t
    end

    o.on("-d","--min-dp NUM",OptionParser::DecimalInteger,"MinDP value to pass to final vcftools recode, defualt #{options[:dp]}") do |t|
      options[:dp] = t
    end

    o.on("-o","--output BASE","Base output path & name, default #{options[:output_base]}") do |p|
      options[:output_base] = File.expand_path p
    end

    o.on("-i","--input GVCF",Array,"Input g.vcf file(s) to GenotypeGVCF") do |g|
      options[:inputs] += g.map {|f| File.expand_path f}
    end

    o.on("-v","--verbose","Output additional verbose information") do |c|
      options[:verbose] = true
    end

    o.on("-D","--debug","Debugging mode, does not actually upload") do |c|
      options[:debug] = true
    end

    o.on("-h","--help","Show this help message") do
      puts o
      exit(0)
    end
    case method
    when :help
      return o.to_s
    else
      o.parse!
    end
  end

  return options
end

def work(jobs)
  passed = true
  workers = jobs.size.times.map do |job|
    Thread.new do
      begin
        prefix = %W/qsub -V -sync y -b y -j y -m e -cwd/
        while job = jobs.pop(true)
          thread = if job[:sge][:threads] && job[:sge][:threads] > 1
                     %W/-pe threaded #{job[:sge][:threads]}/
                   else
                     []
                   end
          mem = if job[:sge][:mem]
                  h=job[:sge][:mem]
                  v=(h/2/job[:sge][:threads]).ceil
                  v = 1 if v <= 1
                  f=(v/2).ceil
                  f = 1 if f <= 1
                  %W/-l h_vmem=#{h}G,virtual_free=#{v}G,mem_free=#{f}G/
                else
                  []
                end
          cmd = prefix + thread + mem + job[:sge][:opts].split(/ /)+ %W/-N #{job[:name]}/ + [job[:cmd].join(" ")]
          if job[:debug]
            puts cmd.inspect
          else
            pid = spawn(job[:env], *cmd,STDOUT=>STDERR)
            pid, status = Process.wait2(pid)
            if nil == status
              STDERR.puts "Unable to start qsub"
              passed = false
            elsif 0 != status.exitstatus
              STDERR.puts "Failed to run #{job[:cmd]}: #{$?}"
              passed = false
            else #clean run, burn in hell log file
              (Dir.glob("#{job[:name]}.po[0-9]*") + Dir.glob("#{job[:name]}.o[0-9]*")).each do |l|
                File.delete(l)
              end
            end
          end
        end #while job
      rescue ThreadError
      end
    end #Thread
  end
  workers.map(&:join)
  return passed
end #workers

def extract_gatk_options(file)
  data = YAML::load_file(file)['DEFAULT']
  if nil == data
    STDERR.puts "Config missing DEFAULT item"
    exit(1)
  end
  data = data.first

  return {
    conf_file:file,
    reference:data[:gatk_ref],
    snp:data[:snp_rod],
    qsub:data[:opts][:qsub_opts] || ""
  }
end

def genotype_gvcfs(gatk_opts,opts)
  STEPS_DIRS_FILES[:genotypegvcf][:files] = STEPS_DIRS_FILES[:genotypegvcf][:ext].call(File.basename(opts[:output_base]))
  cmd = %W/gatk -T GenotypeGVCFs
           -R #{gatk_opts[:reference]}
           -D #{gatk_opts[:snp]}/
  cmd += %W/-nt #{opts[:threads]}/ if opts[:threads] && opts[:threads] > 1
  cmd += opts[:inputs].map {|i| ["-V",i]}.flatten
  cmd += ['-o', STEPS_DIRS_FILES[:genotypegvcf][:files].first ]
  env = {
    "JAVA_MEM_OPTS" => "-Xmx#{opts[:mem]-6}G"
  }
  sge = {
    threads:opts[:threads],
    mem:opts[:mem],
    opts:gatk_opts[:qsub]
  }
  jobs = Queue.new()
  jobs << { name:"genotypeGVCF-#{File.basename(opts[:output_base])}", env:env, sge:sge, cmd:cmd, debug:opts[:debug] }
  passed = false
  Dir.mkdir(STEPS_DIRS_FILES[:genotypegvcf][:dir])
  Dir.chdir(STEPS_DIRS_FILES[:genotypegvcf][:dir]) do
    passed = work(jobs)
  end
  return passed
end

def split_snp_indels(gatk_opts, opts)
  input_file = STEPS_DIRS_FILES[:genotypegvcf][:files].first
  STEPS_DIRS_FILES[:split][:files] = STEPS_DIRS_FILES[:split][:ext].call(File.basename(input_file,".vcf.gz"))
  jobs = Queue.new()
  mem = (opts[:mem]/2).ceil
  STEPS_DIRS_FILES[:split][:files].each_with_index do |f,i|
    cmd = %W/gatk -T SelectVariants -U LENIENT_VCF_PROCESSING
             -R #{gatk_opts[:reference]} -V/
    cmd << File.join("..",STEPS_DIRS_FILES[:genotypegvcf][:dir],input_file)
    cmd += ['-o', f]
    cmd += ['-selectType', STEPS_DIRS_FILES[:split][:types][i]]
    env = {
      "JAVA_MEM_OPTS" => "-Xmx#{mem-6}G"
    }
    sge = {
      threads:1,
      mem:mem,
      opts:gatk_opts[:qsub]
    }
    jobs << { name:"split#{f}-#{File.basename(opts[:output_base])}", env:env, sge:sge, cmd:cmd, debug:opts[:debug] }
  end
  Dir.mkdir(STEPS_DIRS_FILES[:split][:dir])
  passed = false
  Dir.chdir(STEPS_DIRS_FILES[:split][:dir]) do
    passed = work(jobs)
  end
  return passed
end

def vqsr(gatk_opts, opts)
  jobs = Queue.new()
  passed = false
  Dir.chdir(STEPS_DIRS_FILES[:vqsr][:dir]) do
    STEPS_DIRS_FILES[:split][:files].each_with_index do |f,i|
      type = STEPS_DIRS_FILES[:vqsr][:types][i]
      dir = "#{type}s_vqsr"
      Dir.mkdir dir
      STEPS_DIRS_FILES[:vqsr][:files] << File.join(dir,STEPS_DIRS_FILES[:vqsr][:ext].call(File.basename(f,".vcf.gz")))
      cmd = %W/vqsr_vcf.rb -c #{opts[:conf_file]} -i #{f} -m #{type} -o #{dir}/
      sge = {
        threads:1,
        mem:opts[:mem]/2,
        opts:gatk_opts[:qsub]
      }
      jobs << { name:"vqsr#{f}-#{File.basename(opts[:output_base])}", env:{}, sge:sge, cmd:cmd, debug:opts[:debug] }
    end
    passed = work(jobs)
    if passed
      STEPS_DIRS_FILES[:vqsr][:types].each do |type|
        dir = File.join("#{type}s_vqsr","#{Time.now.strftime("%Y%m%d")}_vqsr_vcf_work.")
        Dir.glob("#{dir}[0-9]*").each do |d|
          if File.directory?(d)
            FileUtils.remove_entry_secure(d)
          end
        end
      end
    end
  end # split dir
  return passed
end

def merge_snp_indels(gatk_opts,opts)
  STEPS_DIRS_FILES[:merge][:files] = STEPS_DIRS_FILES[:merge][:ext].call(File.basename(opts[:output_base]))
  cmd = %W/gatk -T CombineVariants -genotypeMergeOptions UNSORTED --assumeIdenticalSamples
           -R #{gatk_opts[:reference]}/
  cmd += %W/-nt #{opts[:threads]}/ if opts[:threads] && opts[:threads] > 1
  cmd += STEPS_DIRS_FILES[:vqsr][:files].map {|i| ["-V",File.join("..",STEPS_DIRS_FILES[:vqsr][:dir],i)]}.flatten
  cmd += ['-o', STEPS_DIRS_FILES[:merge][:files].first ]
  env = {
    "JAVA_MEM_OPTS" => "-Xmx#{opts[:mem]-6}G"
  }
  sge = {
    threads:opts[:threads],
    mem:opts[:mem],
    opts:gatk_opts[:qsub]
  }
  jobs = Queue.new()
  jobs << { name:"merge-#{File.basename(opts[:output_base])}", env:env, sge:sge, cmd:cmd, debug:opts[:debug] }
  passed = false
  Dir.mkdir(STEPS_DIRS_FILES[:merge][:dir])
  Dir.chdir(STEPS_DIRS_FILES[:merge][:dir]) do
    passed = work(jobs)
  end
  return passed
end

def recode(gatk_opts,opts)
  STEPS_DIRS_FILES[:recode][:files] = STEPS_DIRS_FILES[:recode][:ext].call(File.basename(STEPS_DIRS_FILES[:merge][:files].first,".vcf.gz"),opts[:dp])
  input = File.join("..",STEPS_DIRS_FILES[:merge][:dir],STEPS_DIRS_FILES[:merge][:files].first)
  cmd = %W/vcftools --gzvcf #{input} --minDP #{opts[:dp]} --recode-INFO-all
           --recode --keep-INFO-all --stdout/
  if opts[:keep]
    cmd += %W/--keep #{opts[:keep]}/
  end
  cmd += ['>', STEPS_DIRS_FILES[:recode][:files].first ]
  sge = {
    threads:1,
    mem:opts[:mem],
    opts:gatk_opts[:qsub]
  }
  jobs = Queue.new()
  jobs << { name:"recode-#{File.basename(opts[:output_base])}", env:{}, sge:sge, cmd:cmd, debug:opts[:debug] }
  passed = false
  Dir.mkdir(STEPS_DIRS_FILES[:recode][:dir])
  Dir.chdir(STEPS_DIRS_FILES[:recode][:dir]) do
    passed = work(jobs)
  end
  return passed
end

def tabix(gatk_opts,opts)
  STEPS_DIRS_FILES[:tabix][:files] = STEPS_DIRS_FILES[:tabix][:ext].call(File.basename(STEPS_DIRS_FILES[:recode][:files].first,".vcf"))
  cmd = %W/bgzip #{STEPS_DIRS_FILES[:recode][:files].first}/
  sge = {
    threads:1,
    mem:opts[:mem]/2,
    opts:gatk_opts[:qsub]
  }
  jobs = Queue.new()
  jobs << { name:"bgzip-#{File.basename(opts[:output_base])}", env:{}, sge:sge, cmd:cmd, debug:opts[:debug] }
  passed = false
  Dir.chdir(STEPS_DIRS_FILES[:recode][:dir]) do
    passed = work(jobs)
  end
  return passed unless passed
  cmd = %W/tabix -f -p vcf #{STEPS_DIRS_FILES[:tabix][:files].first}/
  jobs << { name:"tabix-#{File.basename(opts[:output_base])}", env:{}, sge:sge, cmd:cmd, debug:opts[:debug] }
  Dir.chdir(STEPS_DIRS_FILES[:recode][:dir]) do
    passed = work(jobs)
  end
  return passed
end

def cleanup()
  STEPS_DIRS_FILES.each do |step,dirs_files|
    next if :tabix == step
    dirs_files[:files].each do |f|
      file = File.join(dirs_files[:dir],f)
      ["",".tbi",".idx"].each do |ext|
        File.delete(file+ext) if File.exists?(file+ext)
      end
      begin
        Dir.rmdir(dirs_files[:dir])
      rescue SystemCallError
      end
    end
  end
end

def main
  options = parsed_opts()
  %w/conf_file/.each do |k|
    unless options[k.to_sym]
      STDERR.puts "Missing #{k} parameter"
      STDERR.puts parsed_opts(:help)
      exit(1)
    end
  end

  gatk_opts = extract_gatk_options(options[:conf_file])

  %w/gatk qsub vcftools bgzip tabix/.each do |b|
    unless which(b)
      STDERR.puts "Unable to find #{b} in PATH"
      exit(1)
    end
  end

  passed = false
  Dir.chdir File.dirname options[:output_base] do
    passed = genotype_gvcfs(gatk_opts,options) &&
      split_snp_indels(gatk_opts, options) &&
      vqsr(gatk_opts, options) &&
      merge_snp_indels(gatk_opts,options) &&
      recode(gatk_opts,options) &&
      tabix(gatk_opts,options) &&
      cleanup()
  end
  unless passed
    exit(1)
  end
end

main()
