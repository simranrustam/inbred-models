#!/bin/bash

set -x

PROJ=$(pwd | rev | cut -d '/' -f-2 | tr '/' '\n' | rev | grep RP)
TSV=$PROJ.tsv

SRP=$(echo $TSV | cut -d '.' -f1)

for SRX in $(cut -f2 $TSV | sed 1d | sort -u )  ; do

  #####################################################
  echo Starting $SRX. Downloading data from SRA and converting to fastq.gz
  #####################################################

  rm ${SRX}*fastq

  for SRR in $(grep -w $SRX $TSV | cut -f3 | sort -u) ; do

    prefetch -a "/usr/local/bin/ascp|/home/mdz/.ascp/aspera-license" $SRR \
    || ( echo $SRR failed download with prefetch ) && \

    mv ~/ncbi/public/sra/${SRR}.sra .

    fastq-dump --split-files ${SRR}.sra && rm ${SRR}.sra

    CNT=$(ls ${SRR}*fastq | wc -l)

    if [ $CNT -eq 1 ] ; then

      cat $SRR.fastq >> $SRX.fastq && rm $SRR.fastq

    fi

    if [ $CNT -eq 2 ] ; then

      cat ${SRR}_1.fastq >> ${SRX}_1.fastq && rm ${SRR}_1.fastq

      cat ${SRR}_2.fastq >> ${SRX}_2.fastq && rm ${SRR}_2.fastq

    fi

  done

  pigz -f ${SRX}*.fastq

  #####################################################
  # need to map with BWA MEM before doing the GRIDSS pipeline
  # skewer trimming
  #####################################################

  if [ $CNT -eq 2 ] ; then

    skewer -q 20 -t 8 ${SRX}_1.fastq.gz ${SRX}_2.fastq.gz

    bwa mem -t 8 \
    /mnt/md0/inbred_models/ref/Mus_musculus.GRCm38.dna_sm.toplevel.fa \
    ${SRX}_1.fastq-trimmed-pair1.fastq ${SRX}_1.fastq-trimmed-pair2.fastq > ${SRX}.sam
    samtools sort -@8 -o ${SRX}.bam ${SRX}.sam \
    && rm ${SRX}.sam ${SRX}_1.fastq-trimmed-pair1.fastq ${SRX}_1.fastq-trimmed-pair2.fastq

  else

    skewer -q 20 -t 8 ${SRX}_1.fastq.gz

    bwa mem -t 8 \
    /mnt/md0/inbred_models/ref/Mus_musculus.GRCm38.dna_sm.toplevel.fa.gz \
    ${SRX}.fastq-trimmed.fastq > ${SRX}.sam
    samtools sort -@8 -o ${SRX}.bam ${SRX}.sam \
    && rm ${SRX}.sam ${SRX}.fastq-trimmed.fastq

  fi

  ../sw/sambamba-0.7.0-linux-static markdup -p -t 8 ${SRX}.bam ${SRX}.nodup.bam \
  && mv ${SRX}.nodup.bam ${SRX}.bam \
  && mv ${SRX}.nodup.bam.bai ${SRX}.bam.bai

done

ls *bam | parallel samtools index {}

featureCounts --ignoreDup -Q 30 -T 8 -F SAF -a /mnt/md0/inbred_models/ref/1e6.bed.saf -o $SRP.1e6.tsv SRX*.bam
sed 1d ${SRP}.1e6.tsv | cut -f1,7- > ${SRP}.1e6_fmt.tsv

featureCounts --ignoreDup -Q 30 -T 8 -F SAF -a /mnt/md0/inbred_models/ref/1e5.bed.saf -o $SRP.1e5.tsv SRX*.bam
sed 1d ${SRP}.1e5.tsv | cut -f1,7- > ${SRP}.1e5_fmt.tsv

featureCounts --ignoreDup -Q 30 -T 8 -F SAF -a /mnt/md0/inbred_models/ref/1e4.bed.saf -o $SRP.1e4.tsv SRX*.bam
sed 1d ${SRP}.1e4.tsv | cut -f1,7- > ${SRP}.1e4_fmt.tsv

featureCounts --ignoreDup -Q 30 -T 8 -F SAF -a /mnt/md0/inbred_models/ref/1e3.bed.saf -o $SRP.1e3.tsv SRX*.bam
sed 1d ${SRP}.1e3.tsv | cut -f1,7- > ${SRP}.1e3_fmt.tsv

featureCounts --ignoreDup -Q 30 -T 8 -F SAF -a /mnt/md0/inbred_models/ref/Mus_musculus.GRCm38.98.gtf.saf -o $SRP.genes.tsv SRX*.bam
sed 1d ${SRP}.genes.tsv | cut -f1,7- > ${SRP}.genes_fmt.tsv

