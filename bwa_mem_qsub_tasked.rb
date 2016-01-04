#!/usr/bin/env ruby

require 'tmpdir'
require 'tempfile'
require 'uri'

module URI
  class O3 < Generic
    USE_REGISTRY = true
  end
  @@schemes['O3'] = O3
end

args = ARGV.clone
tmp_base = args.shift
output_base = args.shift
reference = args.shift
index = (ENV['SGE_TASK_ID']||1).to_i - 1
threads = (ENV['NSLOTS']) || 1

groups = []

while a = args.shift
  data = {}
  if "paired" == a then
    data[:mode] = "paired"
  else
    data[:mode] = "single"
  end
  data[:tag] = args.shift
  data[:inputs] = ["#{args.shift}"]
  data[:inputs] << "#{args.shift}" if "paired" == data[:mode]

  groups << data
end

data = groups[index]
tag = data[:tag]
output = File.join(output_base,"#{index}.bam")

delete_files = []
data[:inputs].map! do |i|
  if i =~ /^(https?|o3):\/\//
    u = URI(i)
    base = File.basename(u.path)
    f = Tempfile.new(base)
    f.close
    tmp_file = f.path + File.extname(base)
    f.unlink
    delete_files << tmp_file
    cmd = case u.scheme
    when /^https?/
      # curl it
      %W/curl --retry 3 -f -o #{tmp_file} #{i}/
    when /^o3/
      # swift it
      (container,object)= u.path.scan(/^\/([^\/]+)\/(.*)/)[0]
      %W/swift download -R 3 -o #{tmp_file} #{container} #{object}/
    end
    pid = spawn(*cmd,STDOUT=>STDERR)
    pid, status = Process.wait2(pid)
    if nil == status
      raise "Unable to start object download"
    elsif 0 != status.exitstatus
      raise "Failure downloading object: #{$?}"
    end
    tmp_file
  else
    i
  end
end

pre = ""
post = ""
data[:inputs].map! do |i|
  #if i =~ /^http.*/ && i =~ /\.xz$/
    #pre ="(t=`mktemp`;"
    #post = " && rm ${t})"
    #(account,cont,obj)=i.scan(/^https:\/\/o3.omrf.org\/v1\/([^\/]+)\/([^\/]+)\/(.*)/)[0]
    #"<(swift download -R 3 #{cont} #{obj} ; xzcat #{obj}  || rm -f ${t})"
  if i =~ /\.xz$/
    pre ="(t=`mktemp`;"
    post = " && rm ${t})"
    "<(xzcat #{i} || rm -f ${t})"
  else
    i
  end
end

return_val = -1
Dir.mktmpdir(index.to_s,tmp_base) do |tmp_prefix_dir|
  cmd = "#{pre}bwa mem -v 1 -M -t #{threads.to_i-2} -R \"#{tag}\" #{reference} #{data[:inputs].join(" ")}#{post}"
  cmd += "|samtools view -Shu - | samtools sort -O bam -@ 4 -m 4G -T #{tmp_prefix_dir}/#{index} > #{output}"

  puts cmd
  STDOUT.flush
  if system "/bin/bash", "-o", "pipefail", "-o", "errexit", "-c", cmd
    return_val = 0
  else
    return_val = $?.exitstatus
  end
end

exit return_val
