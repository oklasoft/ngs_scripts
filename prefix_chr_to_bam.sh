#!/bin/bash

set -o pipefail
set -u

. /usr/local/Modules/default/init/bash
module purge
module load samtools/0.1.19
module load picard/1.99

BAM="$1"
DIR=`dirname "${BAM}"`
BASE=`basename "${BAM}" .bam`
NEWBAM="${BASE}_chr_prefix.bam"

cd $DIR || exit 1


samtools view -H "${BASE}.bam" | sed -e '/^@SQ/s/SN\:/SN\:chr/' -e '/^@SQ/s/chrMT/chrM/' > "${BASE}_new_header.txt"

if [ "$?" -ne "0" ]; then
  echo -e "Failed to make new header"
  exit 1
fi

picard ReplaceSamHeader INPUT="${BASE}.bam" OUTPUT="${NEWBAM}" HEADER="${BASE}_new_header.txt"

if [ "$?" -ne "0" ]; then
  echo -e "Failed to make new bam"
  exit 1
fi

samtools index "${NEWBAM}"

if [ "$?" -ne "0" ]; then
  echo -e "Failed to make index new bam"
  exit 1
fi

rm "${BASE}_new_header.txt"
