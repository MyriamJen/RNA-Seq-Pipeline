# General RNA-Seq analysis
In the following the basic RNA-Seq analysis pipeline used for this thesis will be shown. The code will show examples but not the original code used 
for the analysis as this would use too much space. All RNA samples were sequenced as paired-end samples. The processing for this basic part was mainly
done on Spartan and the original slurm files are stored on Spartan. 

## 1. Checking the sequence quality using FastQC
First, the overall sequencing quality of the samples was assessed using FastQC. All samples are run in a loop for the initial FastQC analysis. 
For this a file called allsamples.txt is created, where all the files are named without their suffix. Their suffix was added to their variable name and
their output directed to a new folder “1_FastQC” that branches from the parent directory holding the files. The quality assessment was done manually.

```bash
for i in $(cat allsamples.txt);
do fastqc $i".fq.gz" > 1_FastQC/$i;
done
#other options:
# -t *number of threads* (should be max. 6 per 32bit PC)
# --contaminants *own list of possible sequences*
# --adapters *own list of adapters used (not necessary if commercially available kit was used)*
```

## 2. STAR-align to align reads
Subsequently, the reads were aligned to the reference genome version 28 of Plasmodium falciparum 3D7. The reference genome was obtained from the 
PlasmoDB website. For this first, the index file was generated. Then the allsamples.txt file was modified to only show the root of the file names without
their directional information (allsample.txt). The output, samfiles and several log files, were directed into “2_STAR-align” folder. As STAR trims the
files automatically, no extra trimming is necessary. The files were directly turned into sorted bam files, which are indexed using samtools.

```bash
#Generation of the index files
STAR --runThreadN 6 --runMode genomeGenerate --genomeDir 3D7_index_files \
--genomeFastaFiles PlasmoDB-28_Pfalciparum3D7_Genome.fasta \
--sjdbGTFfile PlasmoDB-28_Pfalciparum3D7.gff --sjdbGTFtagExonParentTranscript Parent --genomeSAindexNbases 11 --sjdbOverhang 100 
#Aligning the reads to the indexed genome
for i in $(cat allsample.txt); do
STAR --runThreadN 8 --genomeDir 3D7_index_files/ \
--readFilesIn $i'_1.fq' $i'_2.fq' \
--outFileNamePrefix 2_STAR-align/$i --outSAMtype BAM SortedByCoordinate; 
done
samtools index 2_STAR-align/*.bam
```

## 3. Cufflinks to calculate FPKM values
Next Cufflinks was used to calculate the FPKM values for every gene. The results are stored in folders named per sample in the directory “Cufflinks”. 
The transcripts.gtf output file is the one that contains the desired information in the attributes column (column 8). Additionally, Cuffdiff was run, 
Cuffdiff results were used for information on differential expression between parasite stages. The results are stored in folder “Cuffdiff”.

```bash
for i in $(cat allsample.txt); do
cufflinks -o Cufflinks/$i -p 8 --compatible-hits-norm --library-type fr-secondstrand -G PlasmoDB-28_Pfalciparum3D7.gff 3_Bamfiles/$i'.bam';
done
#include --compatible-hits-norm to normalise the single files like the Cuffdiff normalisation
#Example for Cuffdiff usage on male and female gametocyte data
#depends on the premis that the files differ only in their _M and _F suffix
for i in $(cat name.txt); do
cuffdiff -o Cuffdiff/$i -p 8 PlasmoDB-28_Pfalciparum3D7.gff 3_Bamfiles/$i'_F.bam' 3_Bamfiles/$i'_M.bam';
done
```

## 4. Generate directional Bigwig files for visualisation
In a further step, the bam files are transformed into bigwig files, that can be visualised in IGV or used for further analysis of the coverage data. 
The package used for this purpose is deeptools and the subcommand is bamCoverage. The resulting files are stored in “5_Bigwigfiles”. Again, the 
allsample.txt file is used for the loop.

```bash
#Generate Bigwig files
for i in $(cat allsample.txt); do
bamCoverage -p 8 --normalizeUsing 'RPKM' --minMappingQuality 28 --ignoreDuplicates -b 3_Bamfiles/$i'.bam' -o 5_Bigwig/$i'_rpkm_ndq.bw';
done
```

## 5. Generate a MultiBamSummary
For quality control of the samples, a multibamsummary can be generated. This shows if replicates cluster together and if different samples cluster 
away from each other. This was done for all merged files and finally also for the runs altogether to see if the data is reproducible. The bam files 
that are compared are the input of the function, the output is a npz array, that can be used to plot a heatmap or PCA plot. The output files are all 
stored in the “6_MulitBamSummary” folder. This analysis led to the exclusion of D4 and D11 Male data from Run_2020_macrogen. Afterwards, the feasible 
replicates were fused using samtools merge.

```bash
multiBamSummary bins -p 4 --minFragmentLength 50 --maxFragmentLength 3000 -b 3_Bamfiles/*.bam \
-out 6_MultiBamSummary/*filename*-summary-merge.npz
plotCorrelation -in 6_MultiBamSummary/*filename*-summary-merge.npz -o 6_MultiBamSummary/*filename*-summary-merge-heatmap.pdf \
-c spearman -p heatmap --plotNumbers --plotHeight 42 --plotWidth 29.7 --smartlabels \
--outFileCorMatrix 6_MultiBamSummary/*filename*-summary-merge-heatmap.txt
plotPCA -in 6_MultiBamSummary/*filename*-summary-merge.npz -o 6_MultiBamSummary/*filename*-summary-merge-pca.pdf
```

## 6. Merging the Bamfiles
The bam files of the replicates were merged. The resulting files are then stored in the folder “3_Bamfiles_merged_sorted”. The code is an example 
for merging two bam files. Afterwards, from the merged files again bw files and Cufflinks/Cuffdiff analysis can be performed.

```bash
samtools merge *filename merged*.bam *filename1*.bam *filename2*.bam 
samtools sort *filename merged*.bam > *filename merged*_sorted.bam
samtools index *filename merged*_sorted.bam
