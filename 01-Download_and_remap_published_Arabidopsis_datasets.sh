# This file contains Bash code to download and remap the following publicly available datasets (all obtained from 2 weeks old A.thaliana seedlings of Col-0 ecotype):
# 1) ONT Direct RNA-seq (Parker et al. 2020 - PMID 31931956);
# 2) CAGE-seq (Thieffry et al. 2020 - PMID 32213639);
# 3) PAT-seq (Yu et al. 2019 - PMID 31427469);
# 4) plaNET-seq (Kindgren et al. 2020 - PMID 31863587);
# 5) TSS-seq (Nielsen et al. 2019 - PMID 30707695);
# 6) 3' DRS-seq (Schurch et al. 2014 - PMID 24722185);
# 7) pNET-seq (Zhu et al. 2018 - PMID 30374093);
# 8) chrRNA-seq (Jia et al. 2020 - PMID 32541953);
# 9) Stranded RNA-seq (Kohnen et al. 2016 - PMID 27923878);
# 10) TIF-seq (Thomas et al. 2020 - PMID 32444691);


##### Prerequisites -----------------------------------------

# Download the A.thaliana genome:
genome_dir="." # change to an appropriate directory
wget ftp://ftp.ensemblgenomes.org/pub/plants/release-26/fasta/arabidopsis_thaliana/dna/ -O ${genome_dir}/TAIR10.fa

# Generate genomic index for the STAR aligner:
STAR --runMode genomeGenerate --genomeDir ${genome_dir} --genomeFastaFiles ${genome_dir}/TAIR10.fa --runThreadN 4

# Find coordinates of short non-coding RNA genes in Araport11 (and extend them by 100 bp each side):

wget https://www.arabidopsis.org/download_files/Genes/Araport11_genome_release/Araport11_GFF3_genes_transposons.201606.gff.gz
zcat Araport11_GFF3_genes_transposons.201606.gff.gz | sed '/^##/d' | awk 'BEGIN{OFS="\t"}{if ($3=="rRNA" || $3=="tRNA" || $3=="snRNA" || $3=="snoRNA") print $1, $4-100, $5+100, $3, ".", $7}' | sed 's/^Chr//;s/^M/Mt/;s/^C/Pt/' | sort -k1,1 -k2,2n > ${genome_dir}/rRNA_tRNA_snRNA_snoRNA.bed

############################################################
#####               ONT Direct RNA-seq                 #####
############################################################

# Download raw ONT archives from ENA (accession PRJEB32782):
prefix="ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR376"
wget ${prefix}/ERR3764345/col0_nanopore_drs_1.tar.gz
wget ${prefix}/ERR3764347/col0_nanopore_drs_2a.tar.gz
wget ${prefix}/ERR3764348/col0_nanopore_drs_2b.tar.gz
wget ${prefix}/ERR3764349/col0_nanopore_drs_3.tar.gz
wget ${prefix}/ERR3764351/col0_nanopore_drs_4.tar.gz

# Extract subfolders with FASTQ files:
for file in *tar.gz; do tar xvf $file --wildcards '*fastq/pass*'; done

# Concatenate FASTQ files belonging to the same sample:
find ./*FAH45730* -type f -name "*fastq" -exec cat {} > ont_rep1.fq \;
find ./*FAH77434* -type f -name "*fastq" -exec cat {} > ont_rep2_part1.fq \;
find ./*FAH59362* -type f -name "*fastq" -exec cat {} > ont_rep2_part2.fq \;
cat ont_rep2_part?.fq > ont_rep2.fq && rm ont_rep2_part?.fq
find ./*FAH83697* -type f -name "*fastq" -exec cat {} > ont_rep3.fq \;
find ./*FAH83552* -type f -name "*fastq" -exec cat {} > ont_rep4.fq \;

# Compress the final FASTQ files:
for file in ont_rep{1..4}.fq; do gzip $file; done

# Align long reads using Minimap2 v2.17:
for file in ont_rep{1..4}.fq.gz; do echo $file && minimap2 -t 8 -ax splice -k 14 -L --cs --secondary=no -G 10000 ${genome_dir}/TAIR10.fa $file > ${file/fq.gz/sam}; done

# Convert SAM to sorted BAM:
for file in ont*sam; do echo $file && samtools view -hu $file | samtools sort - -o ${file/sam/bam}; done

# Filter out unmapped reads and reads with low MAPQ:
for file in ont*bam; do echo $file && samtools view -hb -q 10 -F 4 $file > ${file/.bam/_final.bam}; done

# Use the ont*final.bam files as input for TranscriptomeReconstructoR

############################################################
#####                    CAGE-seq                      #####
############################################################

# Download CAGE-seq data from SRA (accession GSE136356):
echo -e "SRR10045003\tcage_rep1\nSRR10045004\tcage_rep2\nSRR10045005\tcage_rep3" > cage_acc.txt

while IFS="\t" read -r line || [[ -n "$line" ]]; do acc="$(cut -f1 <<< $line)" && fn="$(cut -f2 <<< $line)" && echo $acc "->" $fn && fastq-dump --gzip $acc && mv ${acc}.fastq.gz ${fn}.fq.gz; done < cage_acc.txt

# Raw CAGE-seq reads are expected to start with non-template sequence XXX-CAGCAG-G (barcode + EcoP15I site + reverse transcription of the cap);
# Thus, retain only bases 11-35 and also skip low quality reads:
for file in cage*fq.gz; do echo $file && zcat $file | fastx_trimmer -f 10 -l 35 | fastq_quality_filter -q 30 -p 50 -z -o ${file/.fq.gz/_trim.fq.gz}; done

# Align CAGE-seq reads using STAR v2.5.2b:
for file in cage*trim.fq.gz; do echo $file && STAR --genomeDir ${genome_dir} --readFilesIn $file --runThreadN 4 --outFileNamePrefix ${file/trim.fq.gz/} --outSAMmultNmax 1 --alignEndsType Extend5pOfRead1 --readFilesCommand zcat --outSAMtype BAM Unsorted; done

for file in *Aligned*; do mv $file ${file/_Aligned.out/}; done
rm -r *STARtmp *out *tab

# Sort BAM files and remove low MAPQ reads:
for file in cage*bam; do echo $file && samtools view -hq 10 $file | samtools sort - -o ${file/.bam/_final.bam}; done

# Use the cage*final.bam files as input for TranscriptomeReconstructoR

############################################################
#####                     PAT-seq                      #####
############################################################

# Download PAT-seq reads from ENA (accession SRP145554):
prefix="ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR716"
wget ${prefix}/006/SRR7160296/SRR7160296.fastq.gz -O pat_rep1.fq.gz
wget ${prefix}/007/SRR7160297/SRR7160297.fastq.gz -O pat_rep2.fq.gz
wget ${prefix}/009/SRR7160299/SRR7160299.fastq.gz -O pat_rep3.fq.gz

# Raw PAT-seq reads are expected to start with 8bp barcodes followed by oligo-T stretches of variable length:
for file in pat*fq.gz; do echo $file && cutadapt -j 4 -u 8 -q 10 -g "T{150};e=0.05" --max-n 3 $file | cutadapt -a AGATCGGAAGAGC -m 20 -o ${file/.fq.gz/_trim.fq.gz} -; done

# Align PAT-seq reads using STAR v2.5.2b:
for file in pat*trim.fq.gz; do echo $file && STAR --genomeDir ${genome_dir} --readFilesIn $file --runThreadN 4 --outFileNamePrefix ${file/trim.fq.gz/} --outSAMmultNmax 1 --alignEndsType Extend5pOfRead1 --readFilesCommand zcat --outSAMtype BAM Unsorted; done

for file in *Aligned*; do mv $file ${file/_Aligned.out/}; done
rm -r *STARtmp *out *tab

# Sort BAM files and remove low MAPQ reads:
for file in pat*bam; do echo $file && samtools view -hq 10 $file | samtools sort - -o ${file/.bam/_final.bam}; done

# Use the pat*final.bam files as input for TranscriptomeReconstructoR

############################################################
#####             plaNET-seq and pNET-seq              #####
############################################################

# Download plaNET-seq data from SRA (accession GSE131733):
echo -e "SRR9117170\tplanet_rep1_part1\nSRR9117171\tplanet_rep1_part2\nSRR9117172\tplanet_rep2_part1\nSRR9117173\tplanet_rep2_part2" > planet_acc.txt

while IFS="\t" read -r line || [[ -n "$line" ]]; do acc="$(cut -f1 <<< $line)" && fn="$(cut -f2 <<< $line)" && echo $acc "->" $fn && fastq-dump --split-files $acc && mv ${acc}_1.fastq.gz ${fn}_R1.fq.gz && mv ${acc}_2.fastq.gz ${fn}_R2.fq.gz; done < planet_acc.txt

# Merge HiSeq and MiSeq files for the same sample:
fn1="planet_rep1_part?_1.fastq"
cat ${fn1} | gzip > planet_rep1_R1.fq.gz && rm ${fn1}
fn2="planet_rep1_part?_2.fastq"
cat ${fn2} | gzip > planet_rep1_R2.fq.gz && rm ${fn2}
fn3="planet_rep2_part?_1.fastq"
cat ${fn3} | gzip > planet_rep2_R1.fq.gz && rm ${fn3}
fn4="planet_rep2_part?_2.fastq"
cat ${fn4} | gzip > planet_rep2_R2.fq.gz && rm ${fn4}

# Download pNET-seq data from SRA (accession GSE109974):
echo -e "SRR6661081\tpnet_s1_total\nSRR6661082\tpnet_s2_total\nSRR6661083\tpnet_s3_unph\nSRR6661084\tpnet_s4_unph\nSRR6661085\tpnet_s5_Ser2P\nSRR6661086\tpnet_s6_Ser2P\nSRR6661087\tpnet_s7_Ser5P\nSRR6661088\tpnet_s8_Ser5P" > pnet_acc.txt

while IFS="\t" read -r line || [[ -n "$line" ]]; do acc="$(cut -f1 <<< $line)" && fn="$(cut -f2 <<< $line)" && echo $acc "->" $fn && fastq-dump --gzip --split-files $acc && mv ${acc}_1.fastq.gz ${fn}_R1.fq.gz && mv ${acc}_2.fastq.gz ${fn}_R2.fq.gz; done < pnet_acc.txt

# Process UMIs in PE mode using UMI-Tools v1.0.1:
for f1 in planet*R1.fq.gz pnet*R1.fq.gz; do f2=${f1/_R1/_R2} && echo $f1 $f2 && umi_tools extract --stdin=${f1} --read2-in=${f2} --bc-pattern=NNNN --bc-pattern2=NNNN --stdout=${f1/.fq.gz/_UMI.fq.gz} --read2-out=${f2/.fq.gz/_UMI.fq.gz}; done

# Trim reads using Trim Galore v0.4.3:
for f1 in planet*R1_UMI.fq.gz pnet*R1_UMI.fq.gz; do f2=${f1/_R1/_R2} && echo $f1 $f2 && trim_galore -q 10 -a TGGAATTCTCGG -a2 GATCGTCGGACT --three_prime_clip_R1 4 --three_prime_clip_R2 4 --paired --gzip --length 15 --no_report_file $f1 $f2; done
for file in *val_1*gz; do mv $file ${file/val_1/trim}; done
for file in *val_2*gz; do mv $file ${file/val_2/trim}; done

# Align to TAIR10 in PE mode using STAR v2.5.2b:
for f1 in planet*R1_UMI_trim.fq.gz pnet*R1_UMI_trim.fq.gz; do f2=${f1/_R1/_R2} && echo $f1 $f2 && STAR --genomeDir ${genome_dir} --readFilesIn $f1 $f2 --runThreadN 4 --outFileNamePrefix ${f1/R1_UMI_trim.fq.gz/} --outSAMmultNmax 1 --alignEndsType Extend5pOfReads12 --readFilesCommand zcat --outSAMtype BAM Unsorted; done

for file in *Aligned*; do mv $file ${file/_Aligned.out/}; done
rm -r *STARtmp *out *tab

# Sort BAM files and skip reads not in proper pairs:
for file in planet*bam pnet*bam; do echo $file && samtools view -hb -f 2 $file | samtools sort - -o ${file/.bam/_sorted.bam}; done

# Filter for reads with high MAPQ:
for file in planet*sorted.bam pnet*sorted.bam; do echo $file && samtools view -q 10 $file -o ${file/.bam/_mapq.bam}; done

# Deduplicate (UMI-Tools):
for file in planet*mapq.bam pnet*mapq.bam; do echo $file && samtools index $file && umi_tools dedup --stdin=${file} --stdout=${file/_sorted_mapq.bam/_final.bam} --paired; done

# Use the planet*final.bam files as input for TranscriptomeReconstructoR

# Also convert first bases of plaNET-seq R2 reads to Bedgraph track (for Fig. 3C):
file="temp_planet_R2.bam"; samtools merge - planet_rep?_final.bam | samtools view -h -f 128 > $file

bedtools genomecov -ibam $file -bg -5 -strand - > ${file/R2.bam/fw.bg}

bedtools genomecov -ibam $file -bg -5 -strand + | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,"-"$4}' > ${file/R2.bam/rev.bg}

cat ${file/R2.bam/fw.bg} ${file/R2.bam/rev.bg} | sort -k1,1 -k2,2n | sed '1i track type=bedGraph color=0,100,200 altColor=200,100,0' | gzip > planet_first_base.bedgraph.gz && rm temp*

# Convert whole pNET-seq inserts to Bedgraph track (for Figures 4C, S2 and S3):

file="temp_pnet.bam"; samtools merge $file pnet_rep?_final.bam

samtools view -F 16 -f 64 $file | awk '{print $3"\t"$4"\t"$4+$9}' | sort -k1,1 -k2,2n | bedtools genomecov -i stdin -bg > temp_pnet_fw.bed && samtools view -F 16 -f 128 $file | awk '{print $3"\t"$4"\t"$4+$9}' | sort -k1,1 -k2,2n | bedtools genomecov -i stdin -bg | awk '{print $1"\t"$2"\t"$3"\t"-$4}' > temp_pnet_rev.bed && cat temp_pnet_fw.bed temp_pnet_rev.bed | sort -k1,1 -k2,2n > pnet_whole_insert.bedgraph.gz && rm temp*

############################################################
#####                     TSS-seq                      #####
############################################################

# Download TSS-seq data from SRA (accession GSE113677):
echo -e "SRR7064004\ttss_rep1_part1\nSRR7064005\ttss_rep1_part2\nSRR7064006\ttss_rep1_part3\nSRR7064007\ttss_rep1_part4\nSRR7064008\ttss_rep2_part1\nSRR7064009\ttss_rep2_part2\nSRR7064010\ttss_rep2_part3\nSRR7064011\ttss_rep2_part4" > tss_acc.txt

while IFS="\t" read -r line || [[ -n "$line" ]]; do acc="$(cut -f1 <<< $line)" && fn="$(cut -f2 <<< $line)" && echo $acc "->" $fn && fastq-dump $acc && mv ${acc}.fastq ${fn}.fq; done < tss_acc.txt

# Merge files corresponding to the same sample:
fn1="tss_rep1_part?.fastq"
cat ${fn1} | gzip > tss_rep1.fq.gz && rm ${fn1}
fn2="tss_rep2_part?.fastq"
cat ${fn2} | gzip > tss_rep2.fq.gz && rm ${fn2}

# Trim custom adapter sequence by Trim Galore v0.4.3:
for file in tss*fq.gz; do echo $file && trim_galore --adapter "ATCTCGTATGCCG" $file; done

# Trim UMIs (first 8 nt of each read) by UMI-Tools v1.0.1:
for file in tss*trimmed.fq.gz; do echo $file && umi_tools extract --stdin=${file} --bc-pattern=NNNNNNNN --stdout=${file/.fq.gz/_UMI.fq.gz}; done

# Align TSS-seq reads using STAR v2.5.2b:
for file in tss*UMI.fq.gz; do echo $file && STAR --genomeDir ${genome_dir} --readFilesIn $file --runThreadN 4 --outFileNamePrefix ${file/.fq.gz/_} --outSAMmultNmax 1 --alignEndsType EndToEnd --readFilesCommand zcat; done

for file in *Aligned*; do mv $file ${file/_Aligned.out/}; done
rm -r *STARtmp *out *tab

# Sort SAM files and convert to BAM:
for file in tss*sam; do echo $file && samtools view -hu $file | samtools sort - -o ${file/.sam/_sorted.bam} && rm $file; done

# Filter out rRNA, tRNA and sn/snoRNA reads:
for file in tss*sorted.bam; do echo $file && bedtools intersect -v -abam $file -b ${genome_dir}/rRNA_tRNA_snRNA_snoRNA.bed > ${file/.bam/_filt.bam}; done

# Filter out multimapper reads:
for file in tss*filt.bam; do echo $file && samtools view -h -q 10 $file -o ${file/.bam/_mapq.bam}; done

# Deduplicate on UMIs:
for file in tss*mapq10.bam; do echo $file && umi_tools dedup -I $file -S ${file/trimmed_UMI_sorted_filt_mapq.bam/final.bam}; done

# Merge the replicates:
samtools merge tss_merged_final.bam tss_rep?_final.bam

# Make stranded Bedgraph files (use only the first base of each read):
for str in "+" "-"; do echo $str; [ "$str" = "+" ] && n="fw" || n="rev"; for file in tss*final.bam; do echo $file && bedtools genomecov -ibam $file -bg -5 -strand $str | sort -k 1,1 -k 2,2n > ${file/final.bam/}${n}.bg; done; done

# Merge forward and reverse Bedgraph files corresponding to the same sample (skip singleton positions):
for f1 in tss*fw.bg; do f2=${f1/fw/rev} && out=${f1/fw.bg/fw_rev.bedgraph.gz} && echo $f1 "+" $f2 "=" $out && awk 'BEGIN{OFS="\t"}{print $1,$2,$3,"-"$4}' $f2 | cat $f1 - | awk '{if (sqrt($4*$4)>1) print}' | sort -k1,1 -k2,2n | sed '1i track type=bedGraph color=0,100,200 altColor=200,100,0' | gzip > $out; done

############################################################
#####                     3' DRS-seq                   #####
############################################################

# Download Helicos FASTQ from DDBJ (accession ERP003245):
prefix="ftp://ftp.ddbj.nig.ac.jp/ddbj_database/dra/fastq/ERA223/ERA223202"
wget ${prefix}/ERX267337/ERR294004.fastq.bz2 drs_rep1.fq.bz2
wget ${prefix}/ERX267338/ERR294005.fastq.bz2 drs_rep2.fq.bz2
wget ${prefix}/ERX267339/ERR294006.fastq.bz2 drs_rep3.fq.bz2

# Remove trailing white spaces after read names and convert to GZ:
for file in drs*fq.bz2; do echo $file && bzcat $file | sed 's/ $//' | gzip > ${file/bz2/gz}; done

# Align 3' DRS-seq reads by STAR v2.5.2b:
for file in drs*fq.gz; do echo $file && STAR --genomeDir ${genome_dir} --readFilesIn $file --runThreadN 4 --outFileNamePrefix ${file/.fq.gz/_} --outSAMmultNmax 1 --alignEndsType Local --readFilesCommand zcat --outSAMtype BAM Unsorted; done

for file in *Aligned*; do mv $file ${file/_Aligned.out/}; done
rm -r *STARtmp *out *tab

# Sort BAM files and remove low MAPQ reads:
for file in drs*bam; do echo $file && samtools view -hu -q 10 $file | samtools sort - -o ${file/.bam/_final.bam}; done

# Merge the replicates:
samtools merge drs_merged_final.bam drs*rep?_final.bam

# Make strand-specific Bedgraph files (use only the first base of each read and switch the strand orientation):
for str in "+" "-"; do [ "$str" = "+" ] && n="rev" || n="fw"; for file in drs*final.bam; do sample=${file/_final.bam/} && echo $n $sample && bedtools genomecov -ibam $file -bg -5 -strand $str | sort -k1,1 -k2,2n > ${sample}_${n}.bg; done; done

# Merge forward and reverse Bedgraph files:
f_str="fw"; r_str="rev"; ext=".bg"; for f1 in drs*${f_str}${ext}; do f2=${f1/${f_str}/${r_str}} && out=${f1/${f_str}${ext}/fw_rev.bedgraph.gz} && echo $f1 "+" $f2 "=" $out && awk 'BEGIN{OFS="\t"}{print $1,$2,$3,"-"$4}' $f2 | cat $f1 - | sort -k1,1 -k2,2n | sed '1i track type=bedGraph color=0,100,200 altColor=200,100,0' | gzip > $out; done


############################################################
#####               Stranded RNA-seq                   #####
############################################################

# Download RNA-seq data from SRA (accession GSE81202):
echo -e "SRR3480142\trna_rep1\nSRR3480144\trna_rep2\n" > rna_acc.txt

while IFS="\t" read -r line || [[ -n "$line" ]]; do acc="$(cut -f1 <<< $line)" && fn="$(cut -f2 <<< $line)" && echo $acc "->" $fn && fastq-dump --gzip $acc && mv ${acc}.fastq.gz ${fn}.fq.gz; done < rna_acc.txt

# Align RNA-seq reads by STAR v2.5.2b:
for file in rna*fq.gz; do echo $file && STAR --genomeDir ${genome_dir} --readFilesIn $file --runThreadN 4 --outFileNamePrefix ${file/.fq.gz/_} --outSAMmultNmax 1 --alignEndsType Extend5pOfRead1 --readFilesCommand zcat --clip3pAdapterSeq AGATCGGAAGAGC --outSAMtype BAM Unsorted; done

for file in *Aligned*; do mv $file ${file/_Aligned.out/}; done
rm -r *STARtmp *out *tab

# Sort BAM files, deduplicate on start coordinates, remove intersections with short ncRNA genes and low MAPQ reads:
for file in rna*bam; do echo $file && samtools view -hu -q 10 | samtools sort - -o - | samtools rmdup -s - - | bedtools intersect -v -abam stdin -b ${genome_dir}/rRNA_tRNA_snRNA_snoRNA.bed > ${file/.bam/_final.bam}; done

# Merge the replicates:
samtools merge rna_merged_final.bam rna_rep?_final.bam

# Make strand-specific Bedgraph files (with strand switch):
for str in "+" "-"; do [ "$str" = "+" ] && n="rev" || n="fw"; for file in rna*final.bam; do sample=${file/_final.bam/} && echo $n $sample && bedtools genomecov -ibam $file -bg -split -strand $str > ${sample}_${n}.bg; done; done

# Merge forward and reverse Bedgraph files:
f_str="fw"; r_str="rev"; ext=".bg"; for f1 in rna*${f_str}${ext}; do f2=${f1/${f_str}/${r_str}} && out=${f1/${f_str}${ext}/fw_rev.bedgraph.gz} && echo $f1 "+" $f2 "=" $out && awk 'BEGIN{OFS="\t"}{print $1,$2,$3,"-"$4}' $f2 | cat $f1 - | sort -k1,1 -k2,2n | sed '1i track type=bedGraph color=0,100,200 altColor=200,100,0' | gzip > $out; done


############################################################
#####                   chrRNA-seq                     #####
############################################################

# Download ONT chrRNA-seq data from SRA (accession PRJNA591665):
echo -e "SRR10538401\tchr_rep1\nSRR10538409\tchr_rep2" > chr_acc.txt

while IFS="\t" read -r line || [[ -n "$line" ]]; do acc="$(cut -f1 <<< $line)" && fn="$(cut -f2 <<< $line)" && echo $acc "->" $fn && fastq-dump --gzip $acc && mv ${acc}.fastq.gz ${fn}.fq.gz; done < rna_acc.txt

# Align long reads by Minimap2 v2.17:
for file in chr*fq.gz; do echo $file && minimap2 -t 8 -ax splice -k 14 -L --cs --secondary=no -G 10000 ${genome_dir}/TAIR10.fa $file > ${file/fastq.gz/sam}; done

# Convert SAM to sorted BAM:
for file in chr*sam; do echo $file && samtools view -hu $file | samtools sort - -o ${file/sam/bam}; done

# Filter out unmapped reads and reads with low MAPQ:
for file in chr*bam; do echo $file && samtools view -hb -q 10 -F 4 $file > ${file/.bam/_final.bam}; done

# Merge the replicates:
samtools merged chr_merged_final.bam chr_rep?_final.bam

# Continue in R to restore the strand information of chrRNA-seq reads (see 02-Postprocess_chrRNAseq.R)


############################################################
#####                     TIF-seq                      #####
############################################################

# Download TIF-seq data from SRA (accession GSE129523):
echo -e "SRR8870029\ttif_rep1\nSRR8870030\ttif_rep2\nSRR8870031\ttif_hen2_rep1\nSRR8870032\ttif_hen2_rep2" > tif_acc.txt

while IFS="\t" read -r line || [[ -n "$line" ]]; do acc="$(cut -f1 <<< $line)" && fn="$(cut -f2 <<< $line)" && echo $acc "->" $fn && fastq-dump --split-files $acc && mv ${acc}_1.fastq.gz ${fn}_R1.fq.gz && mv ${acc}_2.fastq.gz ${fn}_R2.fq.gz; done < tif_acc.txt

# Trim adapters from 3' ends of R1 and R2, trim UMI from 5' end of R2, trim polyT stretches from 5' ends of R1 (up to 4 non-T + at least 8 T), discard R1 reads with internal polyT (at least 12 T):
 
for f1 in tif*R1.fq.gz; do f2=${f1/R1/R2} && name=${f1/R1.fq.gz/} && echo $f1 "+" $f2 && cutadapt -j 4 -a "AGGTGACCGG" -A "AGGTGACCGG" -a "AGATCGGAAG" -A "AGATCGGAAG" --nextseq-trim=20 --match-read-wildcards --minimum-length 28 -o cut1 -p cut2 <(zcat $f1) <(zcat $f2) && umi_tools extract -I cut2 --extract-method=string --bc-pattern=NNNNNNNN --read2-in=cut1 --stdout=umi2 --read2-out=umi1 && rm cut1 cut2 && awk -v name=$name 'BEGIN{OFS="\n"}{h=$0; getline seq; getline qh; getline qseq; if (NR==FNR) {match(seq, /^[^T]{0,4}T{8,}/); end=RSTART+RLENGTH; if (RSTART!=0) {if (length(seq)-end > 20) {seq=substr(seq, end, length(seq)); qseq=substr(qseq, end, length(qseq))} else {arr[FNR]=1; next}}; if (seq~/T{12,}/) {arr[FNR]=1; next}; print h, seq, qh, qseq > name"R1_trimmed.fq"} else {if (arr[FNR]==1) next; print h, seq, qh, qseq > name"R2_trimmed.fq"}}' umi1 umi2 && rm umi1 umi2; done

# Gzip trimmed files:
for file in tif*trimmed.fq; do echo $file && gzip $file; done

# Align TIF-seq reads by STAR v2.5.2b:
for f1 in tif*R1_trimmed.fq.gz; do f2=${f1/R1/R2} && echo $f1 $f2 && STAR --genomeDir ${genome_dir} --readFilesIn $f1 $f2 --runThreadN 4 --outFileNamePrefix ${f1/R1_trimmed.fq.gz/} --outSAMmultNmax 1 --alignEndsType EndToEnd --readFilesCommand zcat --outSAMtype BAM Unsorted --alignIntronMax 10000 --alignMatesGapMax 25000; done

for file in *Aligned*; do mv $file ${file/_Aligned.out/}; done
rm -r *STARtmp *out *tab

# Remove low MAPQ reads and reads not in proper pair. Sort and index BAM files:
for file in tif*bam; do echo $file && samtools view -hq 10 -f 2 $file | samtools sort - -o ${file/.bam/_mapq.bam}; done

for file in tif*mapq.bam; do echo $file && samtools index $file; done

# Deduplicate by UMIs:
for file in tif*mapq.bam; do echo $file && umi_tools dedup -I $file -S ${file/.bam/_final.bam} --method cluster --paired; done

# Merge the replicates:
samtools tif_merged_final.bam tif_rep?_final.bam
samtools tif_hen2_merged_final.bam tif_hen2_rep?_final.bam

# Convert BAM to BED (use the whole insert, flip strand info):
for file in tif*final.bam; do echo $file && samtools view -F 16 -f 64 $file | awk 'BEGIN{OFS="\t"}{print $3, $4, $4+$9, ".", ".", "-"}' > temp_tif_rev.bed && samtools view -F 16 -f 128 $file | awk 'BEGIN{OFS="\t"}{print $3, $4, $4+$9, ".", ".", "+"}' > temp_tif_fw.bed && cat temp_tif_fw.bed temp_tif_rev.bed | sort -k1,1 -k2,2n > ${file/bam/bed} && rm temp*; done

# Collapse identical reads:
for file in tif*bed; do echo $file &&cat $file | awk 'BEGIN{OFS="\t"}{arr[$1":"$2":"$3":"$6]++}END{for (key in arr) {split(key, tx, ":"); print tx[1], tx[2], tx[3], ".", arr[key], tx[4]}}' | sort -k1,1 -k2,2n > temp && mv temp $file; done

