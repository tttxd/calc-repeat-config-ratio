# calc-repeat-config-ratio

A Python script for estimating the relative support of repeat-mediated alternative mitochondrial configurations from long-read BAM files using samtools output.

## Requirements
- Python 3
- samtools

## Usage

```bash
python calc_repeat_config_ratio_generic.py \
  -b example.hifi.bam \
  --config1 path1,path2 \
  --config2 path3,path4 \
  --repeat-len number \
  --left-flank number \
  --right-flank number \
  --anchors 20,50,100 \
  --mapq 30 \
  --samtools samtools \
  --label repeatX \
  --out repeatX_ratio.tsv
