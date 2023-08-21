#######
#Diana Proctor
#This is the series of scripts used to call variants on the Candida auris genomes using freebayes, bcftools, and GATK
#We take the intersection of the three tools to identify high confidence positions


### see: isolates.csv
### we start with 141 isolates downloaded on Sept. 12/13 20222 from IL, MA, FL
#get the samples that didn't download with SAM numbers
#!/bin/sh
module load sratoolkit
fasterq-dump -t /lscratch/$SLURM_JOBID SRR12363157
fasterq-dump -t /lscratch/$SLURM_JOBID SRR12526249
fasterq-dump -t /lscratch/$SLURM_JOBID SRR12526248
fasterq-dump -t /lscratch/$SLURM_JOBID SRR12526246
fasterq-dump -t /lscratch/$SLURM_JOBID SRR12526245
fasterq-dump -t /lscratch/$SLURM_JOBID SRR12526244
fasterq-dump -t /lscratch/$SLURM_JOBID SRR12526243

#get the sam numbers
fasterq-dump -t /lscratch/$SLURM_JOBID SAMN10139189
fasterq-dump -t /lscratch/$SLURM_JOBID SAMN10139194
fasterq-dump -t /lscratch/$SLURM_JOBID SAMN10139195
fasterq-dump -t /lscratch/$SLURM_JOBID SAMN10139226
fasterq-dump -t /lscratch/$SLURM_JOBID SAMN10139227
fasterq-dump -t /lscratch/$SLURM_JOBID SAMN10139260
fasterq-dump -t /lscratch/$SLURM_JOBID SAMN10139351
fasterq-dump -t /lscratch/$SLURM_JOBID SAMN10139373
fasterq-dump -t /lscratch/$SLURM_JOBID SAMN10139374
fasterq-dump -t /lscratch/$SLURM_JOBID SAMN10139498
fasterq-dump -t /lscratch/$SLURM_JOBID SAMN10139499
fasterq-dump -t /lscratch/$SLURM_JOBID SAMN10139500
fasterq-dump -t /lscratch/$SLURM_JOBID SAMN10139501
fasterq-dump -t /lscratch/$SLURM_JOBID SAMN10139511
fasterq-dump -t /lscratch/$SLURM_JOBID SAMN10139512
fasterq-dump -t /lscratch/$SLURM_JOBID SAMN10139513
fasterq-dump -t /lscratch/$SLURM_JOBID SAMN10139526
fasterq-dump -t /lscratch/$SLURM_JOBID SAMN10139527
fasterq-dump -t /lscratch/$SLURM_JOBID SAMN10139534
fasterq-dump -t /lscratch/$SLURM_JOBID SAMN10139535
fasterq-dump -t /lscratch/$SLURM_JOBID SAMN10139536
fasterq-dump -t /lscratch/$SLURM_JOBID SAMN10139537
fasterq-dump -t /lscratch/$SLURM_JOBID SAMN10139538
fasterq-dump -t /lscratch/$SLURM_JOBID SAMN10139539
fasterq-dump -t /lscratch/$SLURM_JOBID SAMN10139540
fasterq-dump -t /lscratch/$SLURM_JOBID SAMN10139541
fasterq-dump -t /lscratch/$SLURM_JOBID SAMN10139542
fasterq-dump -t /lscratch/$SLURM_JOBID SAMN10139543
fasterq-dump -t /lscratch/$SLURM_JOBID SAMN10139544
fasterq-dump -t /lscratch/$SLURM_JOBID SAMN10139545
fasterq-dump -t /lscratch/$SLURM_JOBID SAMN10139546
fasterq-dump -t /lscratch/$SLURM_JOBID SAMN10139547
fasterq-dump -t /lscratch/$SLURM_JOBID SAMN10139548
fasterq-dump -t /lscratch/$SLURM_JOBID SAMN10139551
fasterq-dump -t /lscratch/$SLURM_JOBID SAMN10139552
fasterq-dump -t /lscratch/$SLURM_JOBID SAMN13294115
fasterq-dump -t /lscratch/$SLURM_JOBID SAMN14925561
fasterq-dump -t /lscratch/$SLURM_JOBID SAMN14925562
fasterq-dump -t /lscratch/$SLURM_JOBID SAMN14925563
fasterq-dump -t /lscratch/$SLURM_JOBID SAMN14925564
fasterq-dump -t /lscratch/$SLURM_JOBID SAMN14925565
fasterq-dump -t /lscratch/$SLURM_JOBID SAMN14925566
fasterq-dump -t /lscratch/$SLURM_JOBID SAMN14925567
fasterq-dump -t /lscratch/$SLURM_JOBID SAMN14925568
fasterq-dump -t /lscratch/$SLURM_JOBID SAMN14925569
fasterq-dump -t /lscratch/$SLURM_JOBID SAMN14925570
fasterq-dump -t /lscratch/$SLURM_JOBID SAMN14925571
fasterq-dump -t /lscratch/$SLURM_JOBID SAMN14925572
fasterq-dump -t /lscratch/$SLURM_JOBID SAMN14925573
fasterq-dump -t /lscratch/$SLURM_JOBID SAMN14925574
fasterq-dump -t /lscratch/$SLURM_JOBID SAMN14925575
fasterq-dump -t /lscratch/$SLURM_JOBID SAMN14925576
fasterq-dump -t /lscratch/$SLURM_JOBID SAMN14925577
fasterq-dump -t /lscratch/$SLURM_JOBID SAMN14925579
fasterq-dump -t /lscratch/$SLURM_JOBID SAMN14925580
fasterq-dump -t /lscratch/$SLURM_JOBID SAMN14925581
fasterq-dump -t /lscratch/$SLURM_JOBID SAMN14925582
fasterq-dump -t /lscratch/$SLURM_JOBID SAMN14925583
fasterq-dump -t /lscratch/$SLURM_JOBID SAMN14925584
fasterq-dump -t /lscratch/$SLURM_JOBID SAMN14925585
fasterq-dump -t /lscratch/$SLURM_JOBID SAMN14925586
fasterq-dump -t /lscratch/$SLURM_JOBID SAMN14925587
fasterq-dump -t /lscratch/$SLURM_JOBID SAMN14925588
fasterq-dump -t /lscratch/$SLURM_JOBID SAMN14925589
fasterq-dump -t /lscratch/$SLURM_JOBID SAMN14925590
fasterq-dump -t /lscratch/$SLURM_JOBID SAMN14925591
fasterq-dump -t /lscratch/$SLURM_JOBID SAMN14925592
fasterq-dump -t /lscratch/$SLURM_JOBID SAMN14925593
fasterq-dump -t /lscratch/$SLURM_JOBID SAMN14925594
fasterq-dump -t /lscratch/$SLURM_JOBID SAMN14925595
fasterq-dump -t /lscratch/$SLURM_JOBID SAMN14925596
fasterq-dump -t /lscratch/$SLURM_JOBID SAMN14925597
fasterq-dump -t /lscratch/$SLURM_JOBID SAMN14925598
fasterq-dump -t /lscratch/$SLURM_JOBID SAMN14925599
fasterq-dump -t /lscratch/$SLURM_JOBID SAMN14925600
fasterq-dump -t /lscratch/$SLURM_JOBID SAMN14925601
fasterq-dump -t /lscratch/$SLURM_JOBID SAMN14925602
fasterq-dump -t /lscratch/$SLURM_JOBID SAMN14925603
fasterq-dump -t /lscratch/$SLURM_JOBID SAMN14925604
fasterq-dump -t /lscratch/$SLURM_JOBID SAMN14925605
fasterq-dump -t /lscratch/$SLURM_JOBID SAMN14925606
fasterq-dump -t /lscratch/$SLURM_JOBID SAMN14925607
fasterq-dump -t /lscratch/$SLURM_JOBID SAMN14925608
fasterq-dump -t /lscratch/$SLURM_JOBID SAMN15689540
fasterq-dump -t /lscratch/$SLURM_JOBID SAMN15689542
fasterq-dump -t /lscratch/$SLURM_JOBID SAMN15689543
fasterq-dump -t /lscratch/$SLURM_JOBID SAMN15689544
fasterq-dump -t /lscratch/$SLURM_JOBID SAMN15689545
fasterq-dump -t /lscratch/$SLURM_JOBID SAMN15689546
fasterq-dump -t /lscratch/$SLURM_JOBID SAMN15689547
fasterq-dump -t /lscratch/$SLURM_JOBID SAMN15689549
fasterq-dump -t /lscratch/$SLURM_JOBID SAMN15689550
fasterq-dump -t /lscratch/$SLURM_JOBID SAMN15689552
fasterq-dump -t /lscratch/$SLURM_JOBID SAMN15689553
fasterq-dump -t /lscratch/$SLURM_JOBID SAMN15689554
fasterq-dump -t /lscratch/$SLURM_JOBID SAMN15689555
fasterq-dump -t /lscratch/$SLURM_JOBID SAMN15689556
fasterq-dump -t /lscratch/$SLURM_JOBID SAMN15689557
fasterq-dump -t /lscratch/$SLURM_JOBID SAMN15689558
fasterq-dump -t /lscratch/$SLURM_JOBID SAMN15689559
fasterq-dump -t /lscratch/$SLURM_JOBID SAMN15689560
fasterq-dump -t /lscratch/$SLURM_JOBID SAMN15689561
fasterq-dump -t /lscratch/$SLURM_JOBID SAMN15689562
fasterq-dump -t /lscratch/$SLURM_JOBID SAMN15689563
fasterq-dump -t /lscratch/$SLURM_JOBID SAMN15689564
fasterq-dump -t /lscratch/$SLURM_JOBID SAMN15878075
fasterq-dump -t /lscratch/$SLURM_JOBID SAMN15878076
fasterq-dump -t /lscratch/$SLURM_JOBID SAMN15878077
fasterq-dump -t /lscratch/$SLURM_JOBID SAMN15878078
fasterq-dump -t /lscratch/$SLURM_JOBID SAMN15878079
fasterq-dump -t /lscratch/$SLURM_JOBID SAMN15878080
fasterq-dump -t /lscratch/$SLURM_JOBID SAMN15878081
fasterq-dump -t /lscratch/$SLURM_JOBID SAMN15878082
fasterq-dump -t /lscratch/$SLURM_JOBID SAMN15878084
fasterq-dump -t /lscratch/$SLURM_JOBID SAMN15878093
fasterq-dump -t /lscratch/$SLURM_JOBID SAMN16376889
fasterq-dump -t /lscratch/$SLURM_JOBID SAMN16376890
fasterq-dump -t /lscratch/$SLURM_JOBID SAMN16376891
fasterq-dump -t /lscratch/$SLURM_JOBID SAMN16376892
fasterq-dump -t /lscratch/$SLURM_JOBID SAMN16376893
fasterq-dump -t /lscratch/$SLURM_JOBID SAMN16376894
fasterq-dump -t /lscratch/$SLURM_JOBID SAMN16376895
fasterq-dump -t /lscratch/$SLURM_JOBID SAMN16376896
fasterq-dump -t /lscratch/$SLURM_JOBID SAMN16376897
fasterq-dump -t /lscratch/$SLURM_JOBID SAMN16376898
fasterq-dump -t /lscratch/$SLURM_JOBID SAMN16376899
fasterq-dump -t /lscratch/$SLURM_JOBID SAMN16376900
fasterq-dump -t /lscratch/$SLURM_JOBID SAMN16376901
fasterq-dump -t /lscratch/$SLURM_JOBID SAMN16376902
fasterq-dump -t /lscratch/$SLURM_JOBID SAMN16376903
fasterq-dump -t /lscratch/$SLURM_JOBID SAMN16376904
fasterq-dump -t /lscratch/$SLURM_JOBID SAMN16376906
fasterq-dump -t /lscratch/$SLURM_JOBID SAMN16376907
fasterq-dump -t /lscratch/$SLURM_JOBID SAMN16376912


###################################
#Map reads with bwa to the masked genome 
###################################
bwa index reference/Caur.B11891_adqia.v2.Nano.pilon_shortheader.masked.fasta -p B11891_masked
#see this reference for repeatmasker: https://www.repeatmasker.org/faq.html


##make a bwa.sh script that we will run as a swarm
cat << EOF > bwa.sh
	#!/usr/bin/bash
	READ1=\$1 
	READ2=\$2
	OUTPUT=\$3
	module load bwa
	bwa mem -M -t 16 /data/$USER/cauris_genomes/reference/B11891_masked \${READ1} \${READ2} > alignments/\${OUTPUT}
EOF

bash READ1 READ2 OUTPUT

​

####run picardv5.0.sh
#!/usr/bin/bash
#read this tutorial on population genomics; https://informatics.fas.harvard.edu/whole-genome-resquencing-for-population-genomics-fastq-to-vcf.html
module load picard
module load samtools

sam_files=($(ls  *.sam ))
sorted_bam_names=($(echo "${sam_files[*]}" | sed 's/.sam/_sorted.bam/g'))
derep_bam_names=($(echo "${sorted_bam_names[*]}" | sed 's/_sorted.bam/_derep.bam/g'))
picard_alignment_metrics=($(echo "${sorted_bam_names[*]}" | sed 's/_sorted.bam/_sorted_picard_alignment_metrics.txt/g'))
picard_derep_text=($(echo "${sorted_bam_names[*]}" | sed 's/_sorted.bam/_sorted_derep_metrics.txt/g'))
derep_bam_names_sorted=($(echo "${derep_bam_names[*]}" | sed 's/_derep.bam/_derep_sorted.bam/g'))

# get total subscripts in an array
total=${#sam_files[*]}

#loop over all the files in the input directories 
for (( i=0; i<=$(( $total -1 )); i++ ))
do
        java -Xmx96g -jar $PICARDJARPATH/picard.jar SortSam I=${sam_files[i]} O=${sorted_bam_names[i]} SORT_ORDER=coordinate VALIDATION_STRINGENCY=SILENT CREATE_INDEX=false
        java -Xmx96g -jar $PICARDJARPATH/picard.jar CollectAlignmentSummaryMetrics I=${sorted_bam_names[i]} R=/data/$USER/cauris_genomes/reference/Caur.B11891_adqia.Nano.pilon.fasta METRIC_ACCUMULATION_LEVEL=SAMPLE METRIC_ACCUMULATION_LEVEL=READ_GROUP O=${picard_alignment_metrics[i]}
        java -Xmx96g -jar $PICARDJARPATH/picard.jar MarkDuplicates I=${sorted_bam_names[i]} O=${derep_bam_names[i]} METRICS_FILE=${picard_derep_text[i]} VALIDATION_STRINGENCY=SILENT 
        java -Xmx96g -jar $PICARDJARPATH/picard.jar SortSam I=${derep_bam_names[i]} O=${derep_bam_names_sorted[i]} SORT_ORDER=coordinate VALIDATION_STRINGENCY=SILENT CREATE_INDEX=false  
        rm ${sorted_bam_names[i]}
        rm ${derep_bam_names[i]}
done




####run picard to add read groups
cat << EOF > addRG.sh
#!/usr/bin/bash
module load picard
module load samtools
BAMIN=\$1 
BAMOUT=\$2
RG=\$3
java -Xmx96g -jar \$PICARDJARPATH/picard.jar AddOrReplaceReadGroups I=\${BAMIN} O=\${BAMOUT} RGID=\${RG} RGLB=\${RG} RGPL="none" RGPU="none" RGSM=\${RG}
EOF


#####index the bam files
mv rg_out ..
cd /data/$USER/cauris_genomes/rg_out

bam_files=($(ls  *.bam ))
# get total subscripts in an array
total=${#bam_files[*]}

#loop over all the files in the input directories 
for (( i=0; i<=$(( $total -1 )); i++ ))
do
echo "samtools index ${bam_files[i]}"
done > samtools.swarm
swarm -f samtools.swarm   --job-name index -t 2 -g 98 --time 2:00:00 --logdir index_swarm_out --module samtools



#### index the rerference genome
samtools faidx  Caur_prinseq_good_iYsK.fasta 
bwa index Caur_prinseq_good_iYsK.fasta 
java -Xmx96g -jar $PICARDJARPATH/picard.jar CreateSequenceDictionary -R Caur_prinseq_good_iYsK.fasta 

###### get the bam list
ls -d /data/$USER/cauris_genomes/freebayes_trouble/alignments_fixed/*bam > bam.list



###################################################### call variants with freebayes
#we will specify a quality score minimum of 30 since the default i think is 20 and a prior analysis revealed that this reduces number of variants significantly
###################################################### 

#first, we want to check on the bam files before proceeding
#if you get an error saying that the number of chromosomes don't match between the sample and the fasta it's likely the bam file is empty
#using this script to generate genome regions for parallelization will identify problematic bam files
#this script is available from the freebayes github
#this section on freebayes was written with Thomas Atkins
split_ref_by_bai_datasize.py -r Caur_prinseq_good_iYsK.fasta.fai -L bam.list



# split the reference into regions to use freebayes-parallel
genome_size=$(awk '{sum+=$2;} END{print sum;}' Caur_prinseq_good_iYsK.fasta.fai)
# because generate_regions takes in the size of a region, we give it regions of 
# genome size / core so that we get a number of regions approx. equal to ncpus 
fasta_generate_regions.py Caur_prinseq_good_iYsK.fasta.fai $(($genome_size/36)) > fb_regions.txt



######### run freebayes parallel 
#!/usr/bin/bash
module load freebayes
freebayes-parallel fb_regions.txt \
  36 \
  -f Caur_prinseq_good_iYsK.fasta -p 1 \
  --min-base-quality 30 \
  -L bam.list \
  > freebayes_adqia_140samples_2023-04-28.vcf


###################################################### call variants with bcftoools
#mpileup and then bcftools call
###################################################### 

cat << EOF > bcftools.sh
#!/usr/bin/bash
module load bcftools
#[+] Loading samtools 1.13 
bcftools mpileup -Ov -f Caur_prinseq_good_iYsK.fasta  -b bam.list -B  -Q 30 | bcftools call -mv --ploidy 1 -o bcftools_adqia_shortfasta_140samples_mymask_2023-04-28.vcf
EOF

sbatch --mem=128g --cpus-per-task 32 --time 24:00:00 bcftools.sh

###################################################### call variants with gatk
#create a gatk haplotype caller on individual samples
###################################################### 

mapfile -t bam_names < bam.list
cut -d/ -f 6 bam.list > samples
mapfile -t samples < samples
gatk_vcf_out=($(echo "${samples[*]}" | sed 's/.bam/.gatk_HC.g.vcf/g'))
total=${#bam_names[*]}
for (( i=0; i<=$(( $total -1 )); i++ ))
do	
echo "gatk --java-options '-Xmx8g' HaplotypeCaller -R Caur_prinseq_good_iYsK.fasta -I ${bam_names[i]} -O gatk_out/${gatk_vcf_out[i]} -ploidy 1 -ERC GVCF"
done >gatk.swarm 


swarm -f gatk.swarm --job-name gatk -t 16 -g 96  --time 24:00:00 --module GATK --logdir scripts/gatk_swarm_out



#############################################joint genotyping
#nice reference: https://www.programmersought.com/article/10963095137/

#create a file called: intervals.list containing the following

########################group the vcsf so that we can perform joint genotyping

#create the intervals list
cat /data/$USER/cauris_genomes/Caur_prinseq_good_iYsK.fasta | grep '>' | grep -Eo '[^>]*' | grep -Eo '^[^ ]*' > intervals.list

nano gatk2.sh
#!/usr/bin/bash 
#140 samples
#The Genome Analysis Toolkit (GATK) v4.2.4.1; 
module load GATK
gatk GenomicsDBImport -V *vcf
-L intervals.list \
    --genomicsdb-workspace-path my_database 

#cleanup
mkdir haplotype_caller_out
	mv *gatk_HC*  mkdir haplotype_caller_out
	
############ GenotypeGVCFs
gatk --java-options "-Xmx4g" GenotypeGVCFs \
   -R /data/$USER/cauris_genomes/freebayes_trouble/alignments_fixed/Caur_prinseq_good_iYsK.fasta \
   -V gendb://my_database \
   -O gatk140.GenotypeGVCFs.adqia.shortheader_2023-04-28.vcf.gz
 
 
#################################### 
#normalize variants and filter
#this section on normalization was written with Thomas Atkins
##################################### 
cd /data/$USER/cauris_genomes/freebayes_trouble/variants_out


#Add all possible tags to vcf file
for file in *.vcf; do
  bcftools +fill-tags $file -- -t ALL > tmp
  mv tmp $file
done


#make a directory with the raw vcf files in case we screw something up
mkdir -p raw
for file in *vcf; do
  cp $file $(pwd)/raw/
done

#get the number of SNPs in the raw calls
bcftools view -H bcftools_adqia_shortfasta_140samples_mymask_2023-04-28.vcf  |wc -l #1190
bcftools view -H freebayes_adqia_140samples_2023-04-28.vcf  |wc -l #284954
bcftools view -H gatk140.GenotypeGVCFs.adqia.shortheader_2023-04-28.vcf  |wc -l #1401

## we can see from the number of freebayes snps that we need to normalize variants
# # https://genome.sph.umich.edu/wiki/Variant_Normalization
for vcf in bcftools_adqia_shortfasta_140samples_mymask_2023-04-28.vcf  freebayes_adqia_140samples_2023-04-28.vcf gatk140.GenotypeGVCFs.adqia.shortheader_2023-04-28.vcf; do
    bcftools norm -f /data/$USER/cauris_genomes/freebayes_trouble/Caur_prinseq_good_iYsK.fasta  $vcf -Ov > tmp
    mv tmp $vcf
done
#Lines   total/split/realigned/skipped:  1190/0/197/0
#Lines   total/split/realigned/skipped:  284954/0/25996/0
#Lines   total/split/realigned/skipped:  1401/0/0/0


#get the number of SNPs in the raw calls
bcftools view -H bcftools_adqia_shortfasta_140samples_mymask_2023-04-28.vcf  |wc -l #1190
bcftools view -H freebayes_adqia_140samples_2023-04-28.vcf  |wc -l #284954
bcftools view -H gatk140.GenotypeGVCFs.adqia.shortheader_2023-04-28.vcf  |wc -l #1401


# filter the SNPs based on three filters:
# min. 5 reads per sample
# minimum quality score 1
# only biallelic SNPs

num_samples=$(wc -l < samples.txt)
depth_filter=$((num_samples * 5))

bcftools view -v snps bcftools_adqia_shortfasta_140samples_mymask_2023-04-28.vcf | bcftools filter -i 'MIN(INFO/DP)>'${depth_filter} | bcftools filter -i 'QUAl>30' | bcftools view -m2 -M2 -v snps  > bcftools_filtered_masked.vcf
bcftools view -v snps freebayes_adqia_140samples_2023-04-28.vcf| bcftools filter -i 'MIN(INFO/DP)>'${depth_filter} | bcftools filter -i 'QUAl>30' | bcftools view -m2 -M2 -v snps  > freebayes_filtered_masked.vcf
bcftools view -v snps gatk140.GenotypeGVCFs.adqia.shortheader_2023-04-28.vcf | bcftools filter -i 'MIN(INFO/DP)>'${depth_filter} | bcftools filter -i 'QUAl>30' | bcftools view -m2 -M2 -v snps  > gatk_filtered_masked.vcf

#get the number of SNPs in the raw calls
bcftools view -H bcftools_adqia_shortfasta_140samples_mymask_2023-04-28.vcf  |wc -l #900
bcftools view -H freebayes_adqia_140samples_2023-04-28.vcf  |wc -l #1293
bcftools view -H gatk140.GenotypeGVCFs.adqia.shortheader_2023-04-28.vcf  |wc -l #934


mkdir filtered
mv *_filtered filtered


###########################################
#combine calllers
##########################################
cd /data/$USER/cauris_genomes/freebayes_trouble/variants_out/filtered

# get the sample names and put them into a file so we can sort freebayes
# and gatk in the same order before combining
bcftools query -l bcftools_filtered_masked.vcf  | sort > samples.txt

# sort the samples of all three outputs
bcftools view -S samples.txt bcftools_filtered_masked.vcf > bcftools_filtered_masked_sorted.vcf
bcftools view -S samples.txt freebayes_filtered_masked.vcf > freebayes_filtered_masked_sorted.vcf 
bcftools view -S samples.txt gatk_filtered_masked.vcf > gatk_filtered_masked_sorted.vcf 


# we need to zip and index the files before passing them to bcftools
mkdir sorted
cp *sorted.vcf sorted
cd /data/$USER/cauris_genomes/freebayes_trouble/variants_out/filtered/sorted

for file in *vcf; do
bgzip $file;
done

for file in *gz; do
tabix $file;
done

# take the 3-way intersection of all variant callers
bcftools isec -c all -n=3 -Ov -p /data/$USER/cauris_genomes/freebayes_trouble/3way\
  bcftools_filtered_masked_sorted.vcf.gz  \
  freebayes_filtered_masked_sorted.vcf.gz \
  gatk_filtered_masked_sorted.vcf.gz

# bcftools isec produces a lot of extra files that we
# we're going to move the bcftools file to a 3_way_masked.vcf
cd /data/$USER/cauris_genomes/freebayes_trouble/3way
cp 0000.vcf  3_way_masked.vcf
rm 0001.vcf
rm 0002.vcf 

cd /data/$USER/cauris_genomes/freebayes_trouble/3way
#select variants for the one sample that should be a control
bgzip 3_way_masked.vcf  
tabix  3_way_masked.vcf.gz 

#select variants for the one sample that should be a control: UC_B11891_Env_2016
bcftools view --samples-file samples2keep.txt --min-ac=1 --no-update 3_way_masked.vcf.gz > adqia.3way.bcftools.vcf #6 SNPs


## we can see that these are 
#convert the adqia vcf to a bed file format
sed -e 's/chr//' adqia.3way.bcftools.vcf| awk '{OFS="\t"; if (!/^#/){print $1,$2-1,$2,$4"/"$5,"+"}}'> adqia.bed 
#remove the SNPs that were present in the negative control
bedtools intersect -v -a 3_way_masked.vcf.gz  -b adqia.bed  -wa -header > 3_way_masked_adqia_masked.vcf 
#bcftools view -H 3_way_masked_adqia_masked.vcf |wc -l #765

mkdir final
mv 3_way_masked_adqia_masked.vcf final/

cd /data/$USER/cauris_genomes/freebayes_trouble/3way/final
bgzip *vcf
tabix *gz

#get the ts/tv ratio
bcftools stats 3_way_masked_adqia_masked.vcf.gz > bcftools_stats.txt

###################################################### 
# convert to a phylip object and make a tree
###################################################### 
mkdir /data/$USER/cauris_genomes/freebayes_trouble/phylogeny
cp /data/$USER/cauris_genomes/freebayes_trouble/3way/final/* /data/$USER/cauris_genomes/freebayes_trouble/phylogeny

cd /data/$USER/cauris_genomes/freebayes_trouble/phylogeny

#convert to phylip
python /data/Segrelab/loading_dock/vcf2phylip/vcf2phylip.py --input 3_way_masked_adqia_masked.vcf.gz 

#########now let's generate trees
phy1=3_way_masked_adqia_masked.min4.phy

module load iqtree
#[+] Loading iqtree  2.1.2


####### get trees for the variant caller specific files
iqtree2 -s $phy1  -T AUTO -B 10000



###################################################### 
# Get coverage
###################################################### 
module load bbtools
#[+] Loading samtools 1.14  ... 
#[+] Loading bbtools  38.87 

#### bwa
bam_files=($(ls  *.bam))
covnames=($(echo "${bam_files[*]}" | sed 's/.bam/.bbtools.coverage.txt/g'))
	total=${#bam_files[*]}

	for (( i=0; i<=$(( $total -1 )); i++ ))
	do
  	      echo "bbtools pileup in=${bam_files[i]} out=${covnames[i]} "
	done > bbtools_coverage_bwa.swarm
bash bbtools_coverage_bwa.swarm



