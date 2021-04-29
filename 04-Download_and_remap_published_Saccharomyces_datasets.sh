# This file contains Bash code to download and remap the following publicly available datasets (all obtained from Saccharomyces cerevisiae wild type strain BY4741 (S288C) grown in rich medium at 30 degrees Celsius):
# 1) ONT Direct RNA-seq (Garalde et al., 2018 - PMID 29334379);
# 2) CAGE-seq (Lu et al., 2019 - PMID 31076411);
# 3) 3'READS (Liu et al., 2017 - PMID 28916539);
# 4) NET-seq (Marquardt et al., 2014 - PMID 24949978);
# 5) NET-seq (Topal et al., 2019 - PMID 31558720);
# 6) TSS-seq (Malabat et al., 2015 - PMID 25905671);
# 7) Lexogen QuantSeq 3' mRNA-seq (Schmid et al., 2018 - PMID 30157438);


##### Prerequisites -----------------------------------------

# Download the S.cerevisiae genome:
index_sc="sacCer_STAR" # change to an appropriate directory
wget https://hgdownload.soe.ucsc.edu/goldenPath/sacCer3/bigZips/sacCer3.fa.gz -O sacCer3.fa.gz
gunzip sacCer3.fa.gz

# Generate genomic index for the STAR aligner:
STAR --runMode genomeGenerate --genomeDir ${index_sc} \
  --genomeFastaFiles sacCer3.fa --runThreadN 4

# Download the S.pombe genome:
prefix="https://www.pombase.org/data/genome_sequence_and_features/genome_sequence/Schizosaccharomyces_pombe"

wget ${prefix}_chromosome_I.fa.gz -O Spo2_I.fa.gz
wget ${prefix}_chromosome_II.fa.gz -O Spo2_II.fa.gz
wget ${prefix}_chromosome_III.fa.gz -O Spo2_III.fa.gz
wget ${prefix}_mitochondrial_chromosome.fa.gz -O Spo2_MT.fa.gz

zcat Spo2_I.fa.gz | sed '1d;1i >I' | 
  gzip > temp && mv temp Spo2_I.fa.gz
zcat Spo2_II.fa.gz | sed '1d;1i >II' | 
  gzip > temp && mv temp Spo2_II.fa.gz
zcat Spo2_III.fa.gz | sed '1d;1i >III' | 
  gzip > temp && mv temp Spo2_III.fa.gz
zcat Spo2_MT.fa.gz | sed '1d;1i >MT' | 
  gzip > temp && mv temp Spo2_MT.fa.gz
zcat Spo2_chr[I|II|III|MT].fa.gz | 
  gzip > Spo2.fa.gz && rm Spo2_chr[I|II|III|MT].fa.gz

# Concatenate Spo2 and SacCer3 genomes:
zcat Spo2.fa.gz sacCer3.fa.gz > Spo2_sacCer3.fa

# Generate STAR index for the combined Spo2_sacCer3.fa:
index_sc_spo="sacCer3_Spo2_STAR" # change to an appropriate directory
STAR --runMode genomeGenerate \
  --genomeFastaFiles Spo2_sacCer3.fa \
  --runThreadN 4 --genomeDir ${index_sc_spo}

# Download gene annotation for sacCer3:
wget http://sgd-archive.yeastgenome.org/curation/chromosomal_feature/archive/saccharomyces_cerevisiae.20170114.gff.gz

# Create BED file with coordinates of known rRNA and tRNA genes:
zcat saccharomyces_cerevisiae.201701014.gff.gz | 
sed -n '1,/##FASTA/p' | 
sed '/^#/d;s/^chrmt/chrM/' | 
awk 'BEGIN{OFS="\t"}{if ($3=="rRNA_gene" || $3=="tRNA_gene") \
  print $1,$4-100,$5+100,$3,100,$7}' > \
  SacCer3_rRNA_tRNA_ext100bp.bed

# Add the whole chrM to the BED file:
echo -e "chrM\t1\t85779\tchrM_fw\t100\t+\nchrM\t1\t85779\tchrM_rev\t100\t-" >> SacCer3_rRNA_tRNA_ext100bp.bed


############################################################
#####               ONT Direct RNA-seq                 #####
############################################################

# Download FASTQ files:
echo -e "SRR6059707\trep1
SRR6059706\trep2
SRR6059712\trep3
SRR6059711\trep4
SRR6059710\trep5" > Garalde2018_acc.txt

awk '{cmd="fastq-dump --gzip "$1" && \
  mv "$1".fastq.gz Garalde2018_"$2".fastq.gz "; \
  system(cmd)}' Garalde2018_acc.txt

# Align to SacCer3 genome with Minimap2:
for file in Garalde2018*fastq.gz; do 
  echo $file && 
  minimap2 -t 8 -ax splice -k 14 -L --cs --secondary=no \
    -G 10000 sacCer3.fa $file > ${file/fastq.gz/sam}; 
done

# Convert SAM to sorted BAM:
for file in Garalde2018*sam; do 
  echo $file && 
  samtools view -hu $file | 
  samtools sort - -o ${file/sam/bam}; 
done

# Filter out unmapped reads and reads with low MAPQ:
for file in Garalde2018*bam; do 
  echo $file && 
  samtools view -hb -q 10 -F 4 $file > ${file/.bam/_final.bam}; 
done


############################################################
#####                    CAGE-seq                      #####
############################################################

# Libraries prepared by the nAnT-iCAGE protocol are expected to produce reads with 9 non-template bases at the 5' end (3 nt for sample barcode + 6 nt for UMI);
# However these CAGE-seq samples do not seem to contain 5' UMIs;

# Download FASTQ files:
echo -e "SRR7633070\trep1
SRR7633069\trep2" > Lu2019_acc.txt

awk '{cmd="fastq-dump --gzip "$1" && \
  mv "$1".fastq.gz Lu2019_"$2".fastq.gz"; \
  system(cmd)}' Lu2019_acc.txt

# Align to SacCer3:
for file in Lu2019*fastq.gz; do 
  echo $file && 
  STAR --genomeDir ${index_sc} --readFilesIn $file \
    --runThreadN 4 --outFileNamePrefix ${file/.fastq.gz/_} \
    --clip3pAdapterSeq AGATCGGAAGAGC --outSAMmultNmax 1 \
    --alignEndsType Extend5pOfRead1 --readFilesCommand zcat \
    --outSAMtype BAM Unsorted; 
done

rm *out *tab; rmdir *STARtmp
for file in *Aligned*; do mv $file ${file/_Aligned.out/}; done

# Sort BAM files and remove low MAPQ reads:
for file in Lu2019*bam; do 
  echo $file && 
  samtools view -hq 10 $file | 
  samtools sort - -o ${file/.bam/_final.bam}; 
done


############################################################
#####                     3'READS                      #####
############################################################

# 5' end: UMI (4 nt) + few non-template T's;
# 3' end: Illumina Small RNA adapter;

# Download FASTQ files from SRA:
echo -e "SRR5276077\trep1
SRR5276079\trep2" > Liu2017_acc.txt

awk '{cmd="fastq-dump --gzip "$1" && \
  mv "$1".fastq.gz Liu2017_"$2".fastq.gz"; system(cmd)}' Liu2017_acc.txt

# Trim 3' Illumina Small RNA adapters:
for file in Liu2017*fastq.gz; do 
  echo $file && 
  trim_galore --small_rna --three_prime_clip_R1 4 --gzip $file; 
done

# Process UMI on 5' end (SE mode):
for file in Liu2017*trimmed.fq.gz; do 
  echo $file && 
  umi_tools extract --stdin=${file} --bc-pattern=NNNN \
    --stdout=${file/.fq.gz/_UMI.fq.gz}; 
done

# Align to SacCer3 in Local mode:
for file in Liu2017*fq.gz; do 
  echo $file && 
  STAR --genomeDir ${index_sc} --readFilesIn $file \
    --runThreadN 4 --outFileNamePrefix ${file/.fq.gz/_} \
    --outSAMmultNmax 1 --alignEndsType Local \
    --readFilesCommand zcat --outSAMtype BAM Unsorted; 
done

for file in *Aligned*; do mv $file ${file/_Aligned.out/}; done
rm -r *STARtmp *out *tab

# Sort BAM files and remove low MAPQ reads:
for file in Liu2017*bam; do 
  echo $file && 
  samtools view -hq 10 $file | 
  samtools sort - -o ${file/.bam/_mapq.bam}; 
done

# Deduplicate on UMI:
for file in Liu2017*mapq.bam; do 
  echo $file && 
  samtools index $file && 
  umi_tools dedup --stdin=${file} \
    --stdout=${file/mapq.bam/final.bam}; 
done


############################################################
#####                     NET-seq                      #####
############################################################

# Custom 3' adapter (ATCTCGTATGCCGTCTTCTGCTTG);
# Marquardt et al., 2014 (PMID 24949978): no UMI;
# Topal et al., 2019 (PMID 31558720):
# - 6 nt UMI at 5' end;
# - Spike-in control: S.pombe RNA (~10%);

# Download FASTQ files:
echo -e "SRR1197973\tMarquardt2014_rep1
SRR1197974\tMarquardt2014_rep2
SRR8503038\tTopal2019_rep1
SRR8503039\tTopal2019_rep2" > NET-seq_acc.txt

awk '{cmd="fastq-dump --gzip "$1" && \
  mv "$1".fastq.gz "$2".fastq.gz"; \
  system(cmd)}' NET-seq_acc.txt

# Process 5' UMIs in Topal2019 data:
for file in Topal2019*fastq.gz; do 
  echo $file && 
  umi_tools extract --stdin=${file} --bc-pattern=NNNNNN \
    --stdout=${file/.fastq.gz/_UMI.fq.gz}; 
done

# Align Marquardt2014 files to the SacCer3 genome:
for file in Marquardt2014*fastq.gz; do 
  echo $file && 
  STAR --genomeDir ${index_sc} --readFilesIn $file \
    --runThreadN 4 --outFileNamePrefix ${file/.fastq.gz/_} \
    --outSAMmultNmax 1 --alignEndsType Extend5pOfRead1 \
    --readFilesCommand zcat --clip3pAdapterSeq ATCTCGTATGCCG \
    --outSAMtype BAM Unsorted; 
done

# Align Topal2019 files to the combined SacCer3 + Spo2 genomes:
for file in Topal2019*UMI.fq.gz; do 
  echo $file && 
  STAR --genomeDir ${index_sc_spo} --readFilesIn $file \
    --runThreadN 4 --outFileNamePrefix ${file/UMI.fq.gz/} \
    --outSAMmultNmax 1 --alignEndsType Extend5pOfRead1 \
    --readFilesCommand zcat --clip3pAdapterSeq ATCTCGTATGCCG \
    --outSAMtype BAM Unsorted; 
done

for file in *Aligned*; do mv $file ${file/_Aligned.out/}; done
rmdir *STARtmp; rm *out *tab

# Sort Marquardt2014 BAM files by coordinates and filter by MAPQ values:
for file in Marquardt2014*bam; do 
  echo $file && 
  samtools view -hq 10 $file | 
  samtools sort - -o ${file/.bam/_mapq.bam}; 
done

# Sort Topal2019 BAM by coordinates, filter by MAPQ values and remove Spo2 alignments (contig name do not start with "chr"):
for file in Topal2019*bam; do 
  echo $file && 
  samtools view -hq 10 $file | 
  awk '{if ($3~/chr/ || $2~/chr/) print}' | 
  samtools sort - -o ${file/.bam/_mapq.bam}; 
done

# Remove mitochondrial, rRNA and tRNA reads:
for file in Marquardt2014*mapq.bam Topal2019*mapq.bam; do 
  echo $file && 
  bedtools intersect -v -abam $file \
    -b SacCer3_rRNA_tRNA_ext100bp.bed > \
    ${file/mapq.bam/final.bam}; 
done

# Merge Topal2019 replicates:
samtools merge Topal2019_merged_final.bam Topal2019_rep*final.bam

# Make stranded Bedgraph files from Topal2019 data 
# (use whole reads, switch the strand orientation):
for str in "+" "-"; do 
  [ "$str" = "+" ] && n="rev" || n="fw"; 
  for file in Topal2019*final.bam; do 
    sample=${file/_final.bam/} && 
    echo $n $sample && 
    bedtools genomecov -ibam $file -bg -split\
      -strand $str > ${sample}_${n}.bg; 
  done; 
done

# Merge forward and reverse Bedgraph files:
f_str="fw"; r_str="rev"; ext=".bg"; 
for file1 in Topal2019*${f_str}${ext}; do 
  file2=${file1/${f_str}/${r_str}} && 
  outfile=${file1/${f_str}${ext}/fw_rev.bedgraph.gz} && 
  echo $file1 "+" $file2 "=" $outfile && 
  awk 'BEGIN{OFS="\t"}{print $1,$2,$3,"-"$4}' $file2 | 
  cat $file1 - | sort -k1,1 -k2,2n | 
  sed '1i track type=bedGraph color=0,100,200 altColor=200,100,0' | 
  gzip > $outfile; 
done


############################################################
#####                     TSS-seq                      #####
############################################################

# 5' adapter sequence: NNNNCGCGCGNN(G)

# Download FASTQ files:
echo -e "SRR1708689\trep1
SRR1708697\trep2
SRR1708701\trep3
SRR1708693\trep4
SRR1708705\trep5
SRR1708709\trep6" > Malabat2015_acc.txt

awk '{cmd="fastq-dump --gzip "$1" && \
  mv "$1".fastq.gz Malabat2015_"$2".fastq.gz"; \
  system(cmd)}' Malabat2015_acc.txt

# Process UMIs:
for file in Malabat2015*fastq.gz; do 
  echo $file && 
  umi_tools extract --stdin=${file} --bc-pattern=NNNNXXXXXXNNX \
    --stdout=${file/.fastq.gz/_UMI.fq.gz}; 
done

# Align to SacCer3 with STAR:
for file in Malabat2015*UMI.fq.gz; do 
  echo $file && 
  STAR --genomeDir ${index_sc} --readFilesIn $file \
    --runThreadN 4 --outFileNamePrefix ${file/UMI.fq.gz/} \
    --outSAMmultNmax 1 --alignEndsType Extend5pOfRead1 \
    --readFilesCommand zcat --clip3pAdapterSeq AGATCGGAAGAG \
    --clip5pNbases 7 --outSAMtype BAM Unsorted; 
done

for file in *Aligned*; do mv $file ${file/_Aligned.out/}; done
rm -r *STARtmp *out *tab

# Sort BAM files and filter by MAPQ values:
for file in Malabat2015*bam; do 
  echo $file && 
  samtools view -hq 10 $file | 
  samtools sort - -o ${file/.bam/_mapq.bam}; 
done

# Deduplicate on UMIs:
for file in Malabat2015*mapq.bam; do 
  echo $file && 
  samtools index $file && 
  umi_tools dedup --stdin=${file} \
    --stdout=${file/.bam/_dedup.bam}; 
done

# Merge replicates:
samtools merge Malabat2015_merged_mapq_dedup.bam Malabat2015_rep*dedup.bam

# Make stranded Bedgraph files (without strand switch):
for str in "+" "-"; do 
  [ "$str" = "-" ] && 
  n="rev" || n="fw"; 
  for file in Malabat*dedup.bam; do 
    sample=${file/_mapq_dedup.bam/} && 
    echo $n $sample && 
    bedtools genomecov -ibam $file -bg -5 -strand $str | 
    sort -k1,1 -k2,2n > ${sample}_${n}.bg; 
  done; 
done

# Merge forward and reverse Bedgraph files for the same sample:
f_str="fw"; r_str="rev"; ext=".bg"; 
for file1 in *${f_str}${ext}; 
  do file2=${file1/${f_str}/${r_str}} && 
  outfile=${file1/${f_str}${ext}/fw_rev.bedgraph.gz} && 
  echo $file1 "+" $file2 "=" $outfile && 
  awk 'BEGIN{OFS="\t"}{print $1,$2,$3,"-"$4}' $file2 | 
  cat $file1 - | sort -k1,1 -k2,2n | 
  sed '1i track type=bedGraph color=0,100,200 altColor=200,100,0' | 
  gzip > $outfile; 
done


############################################################
#####                   3' mRNA-seq                    #####
############################################################

# Download FASTQ files:
echo -e "SRR6423294\tstr1
SRR6423292\tstr2" > Schmid2018_acc.txt

awk '{cmd="fastq-dump --gzip "$1" && \
  mv "$1".fastq.gz Schmid2018_"$2".fastq.gz"; \
  system(cmd)}' Schmid2018_acc.txt

# Trim Illumina adapters from 3' ends and remnants of poly(A) tails from 5' ends:
for file in Schmid2018*fastq.gz; do 
  echo $file && 
  cutadapt -j 4 -q 10 -g "T{150};e=0.05" --max-n 3 $file | 
  cutadapt -a AGATCGGAAGAGC -m 20 \
    -o ${file/.fastq.gz/_trim.fq.gz} -; 
done

# Align to SacCer3:
for file in Schmid2018*trim.fq.gz; do 
  echo $file && 
  STAR --genomeDir ${index_sc} --readFilesIn $file \
    --runThreadN 4 --outFileNamePrefix ${file/trim.fq.gz/} \
    --outSAMmultNmax 1 --alignEndsType Extend5pOfRead1 \
    --readFilesCommand zcat --clip3pAdapterSeq AGATCGGAAGAG \
    --outSAMtype BAM Unsorted; 
done

for file in *Aligned*; do mv $file ${file/_Aligned.out/}; done
rm -r *STARtmp *out *tab

# Sort BAM files and filter by MAPQ values:
for file in Schmid2018*bam; do 
  echo $file && 
  samtools view -hq 10 $file | 
  samtools sort - -o ${file/.bam/_mapq.bam}; 
done

# Merge replicates:
samtools merge Schmid2018_merged_mapq.bam Schmid2018_str*mapq.bam

# Make Bedgraph files (with strand switch):
for str in "+" "-"; do 
  [ "$str" = "+" ] && 
  n="rev" || n="fw"; 
  for file in *mapq.bam; do 
    sample=${file/_mapq.bam/} && 
    echo $n $sample && 
    bedtools genomecov -ibam $file -bg -5 -strand $str | 
    sort -k1,1 -k2,2n > ${sample}_${n}.bg; 
  done; 
done

# Merge forward and reverse Bedgraph files for the same sample:
f_str="fw"; r_str="rev"; ext=".bg"; 
for file1 in *${f_str}${ext}; 
  do file2=${file1/${f_str}/${r_str}} && 
  outfile=${file1/${f_str}${ext}/fw_rev.bedgraph.gz} && 
  echo $file1 "+" $file2 "=" $outfile && 
  awk 'BEGIN{OFS="\t"}{print $1,$2,$3,"-"$4}' $file2 | 
  cat $file1 - | sort -k1,1 -k2,2n | 
  sed '1i track type=bedGraph color=0,100,200 altColor=200,100,0' | 
  gzip > $outfile; 
done

