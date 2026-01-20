#!/bin/bash
inputDir=$1
VCF=$2
bcftools convert --gvcf2vcf --fasta-ref /home/jslinkspark/work/hg38/Homo_sapiens_assembly38.fasta ${inputDir}/${VCF} -Oz -o ${inputDir}/${VCF%.vcf.gz}.vcf.bgz --threads 15
/home/jslinkspark/work/tools/pharmcat/pharmcat_vcf_preprocessor -G -vcf ${inputDir}/${VCF%.vcf.gz}.vcf.bgz

java -jar /home/jslinkspark/work/tools/pharmcat/pharmcat.jar -vcf ${inputDir}/${VCF%.vcf.gz}.preprocessed.vcf.bgz -reporterCallsOnlyTsv --reporter-extended # --reporter-save-html -reporterJson 