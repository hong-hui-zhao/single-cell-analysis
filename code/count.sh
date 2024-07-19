#!/bin/bash
cellranger count --id=086 \
   --fastqs=/home/zhh/downloads/HRR338087 \
   --sample=LUAD \
   --create-bam=false \
   --transcriptome=/home/zhh/refdata-gex-GRCh38-2024-A \
   --nosecondary
   