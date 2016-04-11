#!/usr/bin/env ruby
# encoding: UTF-8
# Copyright (c) 2016 Stuart Glenn, Oklahoma Medical Research Foundation
# Licensed under 3 Clause BSD
# Wrapper to pull together some fastq and alignment metrics for report

# Things for which we are looking and reporting include
# - Percent undetermined for the demuxed flowcell (bcl2fastq)
# <=10
# - Average %Q30 across all lanes for each sample (bcl2fastq)
sleep .5
# =75
# - Average mean Q score across all lanes for each sample (bcl2fastq)
# =30
#
# - Percent pass filter unique reads (picard)
# PCT_PF_UQ_READS
# =75
# - Percent pass filter unique reads aligned to reference (picard)
# PCT_PF_UQ_READS_ALIGNED
# =90
# - Mean target coverage (picard)
# MEAN_TARGET_COVERAGE
# =100
# - Percent target bases > 10x (picard)
# PCT_TARGET_BASES_10X
# =90
# - Percent target bases > 20x (picard)
# PCT_TARGET_BASES_20X
# =85
# - Percent target bases > 30x (picard)
# PCT_TARGET_BASES_30X
# =80
#
# - VerifyBamID (FREELK1/FREELK0) ratio (VerifyBamID)
# =0.90
# verifyBamID --vcf ${vcf_file.vcf.gz}  --bam ${bam_file.bam} --out ${out_name} --ignoreRG --chip-none --maxDepth 1000 --precise
# read then out_name.selfSM

$:.unshift File.dirname(File.realpath(__FILE__))

require 'yaml'
require 'optparse'
require 'fileutils'
require 'uri'
require 'uri-o3'
require 'net/http'
require 'csv'
require 'yaml'

def metric_class(n,c,perc=false)
  Class.new do
    define_method :check do
      c
    end
    define_method :name do
      n
    end

    attr_reader :val

    define_method :initialize do |val|
      @val = val
    end

    def valid?
      if nil == @val
        return false
      elsif nil == check
        ""
      else
        check.call(@val)
      end
    end

    if perc
      define_method :formated_val do |fmt|
        return "NDA" if nil == @val
        sprintf(fmt,(@val * 100.0).round(2))
      end
    else
      define_method :formated_val do |fmt|
        return "NDA" if nil == @val
        sprintf(fmt,@val.round(2))
      end
    end
  end
end


METRICS = {
  :number_of_runs => metric_class("BCL Number of Runs", lambda {|x| x < 2 && x > 0}),
  :avg_of_lanes => metric_class("BCL Demux Average % of Lanes", lambda {|x| x > 0}),
  :num_lanes => metric_class("BCL Demux Number of Lanes", lambda {|x| x >= 1 && x <= 4}),
  :avg_q30 => metric_class("BCL Demux Average Q30", lambda {|x| x >= 75.0}),
  :avg_unknown => metric_class("BCL Demux Average % Unknown",lambda {|x| x <= 10.0}),
  :avg_mean_q => metric_class("BCL Demux Average Mean Q",lambda {|x| x >= 30.0 }),
  :freelk_ratio => metric_class("VerifyBam FreeLK Ratio",lambda {|x| x >= 0.90}),
  :pct_pf_uq_reads => metric_class("Picard % PF UQ Reads",lambda {|x| x >= 0.75}, true),
  :pct_pf_uq_reads_aligned => metric_class("Picard % PF UQ Reads Aligned",lambda {|x| x >= 0.90}, true),
  :mean_target_coverage => metric_class("Picard Mean Target Coverage",lambda {|x| x >= 100}),
  :pct_target_bases_10x => metric_class("Picard % Target Bases 10x",lambda {|x| x >= 0.90}, true),
  :pct_target_bases_20x => metric_class("Picard % Target Bases 20x",lambda {|x| x >= 0.85}, true),
  :pct_target_bases_30x => metric_class("Picard % Target Bases 30x",lambda {|x| x >= 0.80}, true),
}

@opts = {
  :verbose => false,
  :debug => true,
  :conf => nil,
  :bam => nil,
  :vcf => nil
}

def demux_stats(conf,sample_name)
  c = YAML::load_file(conf)[sample_name]
  runs = c.map{|r| demux_stats_from_o3(r[:inputs].first)}.reject {|r| r.empty?}
  results = {
    :number_of_runs => c.size#runs.size
  }
  return results if runs.empty?
  runs.first.keys.each do |k|
    sum = runs.reduce(0.0) {|a,v| a += v[k]}
    if :num_lanes == k
      results[k] = sum
    else
      results[k] = sum/results[:number_of_runs]
    end
  end
  return results
end

def demux_stats_from_o3(object_path)
  token = ENV['OS_AUTH_TOKEN']
  u = URI.parse(object_path)

  http = Net::HTTP.new(u.host, u.port)
  http.use_ssl=true
  http.verify_mode = OpenSSL::SSL::VERIFY_PEER
  http.open_timeout = 2
  http.read_timeout = 15
  req  = Net::HTTP::Head.new(u.request_uri)
  req["X-Auth-Token"] = token
  resp = http.request(req)
  msg = case resp
        when Net::HTTPSuccess
          r = {}
          %w/avg-q30 avg-unknown avg-mean-q avg-of-lanes num-lanes/.each do |q|
            h_key = "x-object-meta-sequence-run-library-demux-#{q}"
            if resp.key?(h_key)
              r[q.gsub(/-/,'_').to_sym] = resp[h_key].to_f.round(2)
            end
          end
          return r
        when Net::HTTPUnauthorized, Net::HTTPForbidden
          "Unauthorized, check OS_AUTH_TOKEN"
        when Net::HTTPNotFound
          "Object #{object_path} not found"
        when Net::HTTPServerError
          "Server Error: #{resp.value}"
        else
          "Unknown error: #{resp.value}"
        end
  $stderr.puts "O3 Error: #{msg}"
  return {}
end

# FREELK1/FREELK0
def verify_bam_ratio(bam_path,vcf_path)
  sm_tsv = verify_bam_sm_path_for(bam_path,vcf_path)
  return nil if nil == sm_tsv
  d = CSV.read(sm_tsv,{:headers=>true,:col_sep=>"\t",:converters=>:numeric})
  return (d[0]['FREELK1']/d[0]['FREELK0']).round(2)
end

def verify_bam_sm_path_for(bam_path,vcf_path)
  b = File.basename(bam_path,".bam")
  v = File.basename(vcf_path,".vcf.gz")
  out = File.join(File.dirname(bam_path),"verify-#{b}-#{v}")
  final = "#{out}.selfSM"
  if File.exist?(final)
    return final
  end
  cmd = %W/verifyBamID --self --precise --ignoreRG --chip-none --maxDepth 1000 --vcf
          #{vcf_path} --bam #{bam_path} --out #{out}/

  begin
    pid = spawn(*cmd,STDOUT=>'/dev/null',STDERR=>'/dev/null')
  rescue Errno::ENOENT
    $stderr.puts "Missing verifyBamID from $PATH"
    return nil
  end
  pid, status = Process.wait2(pid)
  if nil == status
    $stderr.puts "Unable to start verifyBamID"
    return nil
  elsif 0 != status.exitstatus
    $stderr.puts "Failure verifyBamID #{$?}"
    return nil
  end
  %w/depthSM log/.each do |ext|
    f = "#{out}.#{ext}"
    begin
      File.delete(f) if File.exists?(f)
    rescue
    end
  end
  return final
end

def picard_hs(bam_path,hs,sample_name,conf)
  c = YAML::load_file(conf)
  cmd = %W/picard CollectHsMetrics INPUT=#{bam_path} OUTPUT=#{hs}
           VALIDATION_STRINGENCY=LENIENT REFERENCE_SEQUENCE=#{c['DEFAULT'].first[:gatk_ref]}/
  i_file = c[sample_name].first[:interval_file] || c['DEFAULT'].first[:interval_file]
  if nil == i_file
    $stderr.puts "Picard HS metrics requires an interval/bait file"
    return nil
  end
  cmd += %W/BAIT_INTERVALS=#{i_file} TARGET_INTERVALS=#{i_file}/
  begin
    pid = spawn(*cmd,STDOUT=>'/dev/null',STDERR=>'/dev/null')
  rescue Errno::ENOENT
    $stderr.puts "Missing picard from $PATH"
    return nil
  end
  pid, status = Process.wait2(pid)
  if nil == status
    $stderr.puts "Unable to start picard"
    return nil
  elsif 0 != status.exitstatus
    $stderr.puts "Failure picard #{$?}"
    return nil
  end
  return hs
end

# - Percent pass filter unique reads (picard)
# PCT_PF_UQ_READS
# - Percent pass filter unique reads aligned to reference (picard)
# PCT_PF_UQ_READS_ALIGNED
# - Mean target coverage (picard)
# MEAN_TARGET_COVERAGE
# - Percent target bases > 10x (picard)
# PCT_TARGET_BASES_10X
# - Percent target bases > 20x (picard)
# PCT_TARGET_BASES_20X
# - Percent target bases > 30x (picard)
# PCT_TARGET_BASES_30X
def picard_stats(bam_path,conf)
  b = File.basename(bam_path,".bam")
  hs = File.join(File.dirname(bam_path),"#{b}_hs_metrics.txt")
  unless File.exist?(hs)
    picard_hs(bam_path,hs,b,conf)
    unless File.exist?(hs)
      $stderr.puts "No picard hs metrics: #{hs}"
      return {}
    end
  end
  d = ""
  IO.foreach(hs) do |line|
    if "" != d
      d << line
    end
    next unless line =~ /\ABAIT_SET/
    d << line
  end
  c = CSV.parse(d,{:headers=>true,:col_sep=>"\t",:converters=>:numeric})
  r = {}
  %w/PCT_PF_UQ_READS PCT_PF_UQ_READS_ALIGNED MEAN_TARGET_COVERAGE PCT_TARGET_BASES_10X PCT_TARGET_BASES_20X PCT_TARGET_BASES_30X/.each do |k|
    r[k.downcase.to_sym] = c[0][k].round(2)
  end
  return r
end

op = OptionParser.new do |o|
  o.banner = "Usage: #{File.basename(__FILE__)} -c CONF_FILE -b BAM_FILE -v VCF_FILE"
  o.on("-V","--verbose","Enable verbose output") do
    @opts[:verbose] = true
  end
  o.on("-D","--debug","Enable debug mode") do
    @opts[:debug] = true
  end
  o.on("-h","--help","Show this help message") do
    puts o
    exit(0)
  end
  o.on("-c","--conf YAML","Load YAML for analysis conf, for fastq source") do |y|
    @opts[:conf] = File.expand_path(y)
  end
  o.on("-b","--bam BAM","Specify BAM file to parse") do |b|
    @opts[:bam] = File.expand_path(b)
  end
  o.on("-v","--vcf VCF","Specify VCF file to parse for verifyBamID") do |v|
    @opts[:vcf] = File.expand_path(v)
  end
  o.parse!
end

%w/vcf bam conf/.each do |k|
  unless @opts[k.to_sym] && File.exist?(@opts[k.to_sym])
    $stderr.puts "Missing valid #{k} parameter"
    $stderr.puts op.help()
    exit(1)
  end
end

sample_name = File.basename(@opts[:bam],".bam")

stats = demux_stats(@opts[:conf],sample_name)
stats[:freelk_ratio] = verify_bam_ratio(@opts[:bam],@opts[:vcf])
stats.merge!(picard_stats(@opts[:bam],@opts[:conf]))

puts sample_name
puts "=" * sample_name.length
METRICS.each do |k,mc|
  m = mc.new(stats[k])
  pass = m.valid? ? " PASS ✅" : " FAIL ❌"
  printf("    %28s %7s%s\n",m.name,m.formated_val("%7.2f"),pass)
end