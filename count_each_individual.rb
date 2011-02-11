INTERVAL = "/Volumes/hts_raw/scratch/merge_gaffney_go_20110125_freeze_bam_merge/20110208_wileyg_subset_new_interval/OMRF_EA_sureselect_hg19.interval_list"

DATA.readlines.each do |sample|
  sample.chomp!
  bam_path = "/Volumes/hts_raw/analysis/gaffney_go/20110125_freeze/individuals/all/#{sample}/13_final_bam/#{sample}.bam"
  Dir.mkdir(sample)
  Dir.chdir(sample) do
    cmd = "qsub -V -b y -j y -cwd -m e -N #{sample}_coverage gatk -T DepthOfCoverage -R /Volumes/hts_core/Shared/homo_sapiens_37/chr_fixed/b37.fasta -I #{bam_path} -o #{sample} -L #{INTERVAL}"
    puts cmd
    system(cmd)
    sleep(5)
  end
end

__END__
lgs103077
lgs105288
lgs300285
lgs103614
lgs304174
lgs302003
lgs301263
lgs102919
lgs103873
lgs100440
lgs500080
lgs103201
lgs104491
lgs102841
lgs301424
lgs102008
lgs101688
lgs000373
lgs103181
lgs102274
lgs304448
lgs105407
lgs103702
lgs304801
lgs304526
lgs102018
lgs301950
lgs301842
lgs100813
lgs301088
lgs103050
lgs101292
lgs304473
lgs302804
lgs100626
lgs101220
lgs301936
lgs300140
lgs101201
lgs101917
lgs301832
lgs303491
lgs100292
lgs303453
lgs304476
lgs102850
lgs100664
lgs102796
lgs102724
lgs101395
lgs300860
lgs101291
lgs101609
lgs303570
lgs100403
lgs101886
lgs103161
lgs303605
lgs302734
lgs304614
lgs303627
lgs300997
lgs304082
lgs304595
