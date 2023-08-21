#Author: Diana Proctor
#Start Date: April 18, 2022
#Finish Date: May 5, 2022
#Revision Date: August 19, 2022
#v2.0 differs from v1 in where I'm pulling eukaryotic bins; v1 pulled froom metawrap refined; v2 pulls from initial concoct binning
#v3.0 differs from 2.0 in thee inclusion of 160 additional samples and is cleaned up based on input from biowulf, also megahit integrated


#######Note: this protocol assumes the reads were cleaned using Sean's SOP
#It's describe here: https://segrewiki.nhgri.nih.gov/lab/index.php/File:V2_Metagenomic_Pipeline.docx 
#"Reads were processed by trimming adapters using cutadapt [64], 
#removing low-quality reads using prinseq-lite [65] with option “-min_qual_mean 20” and removing host (mouse; GRCm39) reads."
#as a result we don't begin with the knead data steps outlined in the Nature Methods protocol
########Note: this protocol assumes you have conda installedd
#if not, follow theses steps here: https://hpc.nih.gov/apps/python.html#envs
#and add the file 'mysource' to /home/$USER/bin


######################################################### 
#step 0: rename the fastq files to the format from R1.fq and R2.fq to _1.fastq and _2.fastq since it's required for binning
######################################################### 
####0a. gunzip files
for file in *gz
do
	echo "gunzip $file"
done > gunzip.swarm


#####0b. rename files
#note copy the reads you want into $target
target=/data/$USER/metagenomes_mags/data/00_reads/single
cd  $target

read1=($(ls  *_R1.fq ))
read2=($(ls  *_R2.fq ))

#define outputs
R1_out=($(echo "${read1[*]}" | sed 's/_R1.fq/_1.fastq/g'))
R2_out=($(echo "${read2[*]}" | sed 's/_R2.fq/_2.fastq/g'))


#loop over all the files and rename them
total=${#read1[*]}
for (( i=0; i<=$(( $total -1 )); i++ ))
do
	echo "mv ${read1[i]} ${R1_out[i]}"
	echo "mv ${read2[i]} ${R2_out[i]}"
done > rename.sh
bash rename.sh

mkdir single
mv * single
 
######################################################### 
#step 1: download the databases required for protocol
######################################################### 
sinteractive --cpus-per-task=8 --time 36:00:00  --mem=96g --gres=lscratch:400 

#######step 1a: make a conda environment for ncbi-genome-download and download bacterial refseq
source myconda 
conda activate
#conda create -n ncbi-genome-download
conda activate ncbi-genome-download
#conda install -c bioconda ncbi-genome-download

#######step 1b: make a mash sketch of the bacterial refseq database
#use ncbi-genome-download to download the bacterial refseq database, including only complete assemblies
#download data: May 5, 2022
mkdir /data/$USER/metagenomes_mags/data
mkdir /data/$USER/metagenomes_mags/data/databases
cd /data/$USER/metagenomes_mags/data/databases
	ncbi-genome-download bacteria --formats fasta --section refseq --assembly-levels complete 
	mv refseq/ bacteria_refseq
	cd bacteria_refseq/bacteria
	mkdir fastas
	cp **/*fna.gz  fastas
	gunzip fasta/*gz

	#make a mash sketch of bacterial refseq
	module load mash
	#[+] Loading mash  2.3 

	cd /data/$USER/metagenomes_mags
	mash sketch -o refseq.bacteria.msh /data/$USER/metagenomes_mags/data/databases/bacteria_refseq/bacteria/fastas/*fna
	

#######step 1c: #make a mash sketch of the fungal genbank database
#note I used a mash sketch against refseq and it was insufficient to identify taxa for 5 samples
#use ncbi-genome-download to download the fungal genbank database, including all assemblies
#download data: May 5, 2022
mkdir /data/$USER/metagenomes_mags/data/databases/fungi_genbank
cd /data/$USER/metagenomes_mags/data/databases/fungi_genbank
	ncbi-genome-download fungi --formats fasta --section genbank
		#WARNING: Skipping entry, as it has no ftp directory listed: 'GCA_023212695.1' #Sporisorium sorghi (smut fungi)
		#WARNING: Skipping entry, as it has no ftp directory listed: 'GCA_023184635.1' #Umbilicaria deusta (ascomycetes)
		#WARNING: Skipping entry, as it has no ftp directory listed: 'GCA_023213265.1' #Ustilago tritici (smut fungi)


	cd /data/$USER/metagenomes_mags/data/databases/fungi_genbank/genbank/fungi
	mkdir fastas
	cp **/*fna.gz fastas/
	gunzip fastas/*gz
	ls fastas/*fna |wc -l #11001
	#go to this link and see how many assemblies we should have: https://www.ncbi.nlm.nih.gov/assembly = 11004
	#based on a search of the GCA numbers above we won't download those failed sequences

#######step 1d: grab species specific to your dataset 	
	#grab the fungal species from this outbreak
	cp /data/Segrelab/data/all_genomes/Caur_296.Fg_acrtq.spades.pilon.fasta /data/$USER/metagenomes_mags/data/databases/fungi_genbank/genbank/fungi/fastas #C. parapsilopsis
	cp /data/Segrelab/data/all_genomes/Caur_2117.An_acuzh.spades.pilon.fasta /data/$USER/metagenomes_mags/data/databases/fungi_genbank/genbank/fungi/fastas #C. glabrata
	cp /data/Segrelab/data/all_genomes/Caur_283.Ic_acrty.spades.pilon.fasta /data/$USER/metagenomes_mags/data/databases/fungi_genbank/genbank/fungi/fastas #C. parapsilopsis
	
#######step 1e: grab the SMGC fungal species
	cp /data/Segrelab/data/zoo/SMGCe/*fa /data/$USER/metagenomes_mags/data/databases/fungi_genbank/genbank/fungi/fastas 
	
	
#######step 1f: make a mash sketch of the fungal databases
	module load mash
	#[+] Loading mash  2.3 
	cd /data/$USER/metagenomes_mags/data/databases
	mash sketch -o /data/$USER/metagenomes_mags/data/databases/fungi.genbank.msh  /data/$USER/metagenomes_mags/data/databases/fungi_genbank/genbank/fungi/fastas/* #11011


#######step 1g: we will need taxonomy information, associating these fasta files to species
#go to your web browser and get the taxonomy files associated with the genbank downloads
#mv these to your database directory: /data/$USER/metagenomes_mags/data/databases
	#ftp://ftp.ncbi.nih.gov/genomes/genbank/fungi/assembly_summary.txt
	#mv assembly_summary.txt assembly_summary_genbank_fungi.txt
	#manually add a column titled "Reference2" which contains the filename associated with the ftp filepath (e.g., GCA_000151125.2_ASM15112v2)

	#ftp://ftp.ncbi.nih.gov/genomes/refseq/bacteria/assembly_summary.txt
	#mv assembly_summary.txt assembly_summary_bacteria.txt

#######step 1h: I use this file to merge with the mummer/mash output in later steps for species IDs
#go to this repository and copy and paste the contents into a fileNCBI_taxID_list.txt in /data/$USER/metagenomes_mags/data/databases
	#https://github.com/Mangul-Lab-USC/db.microbiome/blob/18a8228596e373be789bf706f63ac8ca4b99b17f/Fungi/code/NCBI_taxID_list.txt

#######step 1i: download the checkM database, expand the tarball and cleanup
mkdir /data/$USER/metagenomes_mags/data/databases/MY_CHECKM_FOLDER
cd /data/$USER/metagenomes_mags/data/databases/MY_CHECKM_FOLDER
	wget https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz
	tar -xvzf checkm_data_2015_01_16.tar.gz
	rm checkm_data_2015_01_16.tar.gz 

#######step 1j: download a dependency from  Alex Almeida's github
cd /home/$USER/bin
git clone https://github.com/Finn-Lab/MGS-gut.git

	
#we will need to download gatkdb, eukccdb, but we will do this later

#################################################### 
#step 2: spades single sample assemblies
#spades.py --meta -1 {input.fwd} -2 {input.rev} -o
#{params.outdir} --threads {threads} -m {resources.mem}
#launched on May 5, 2022 for 104 samples at 1:27 pm
#80 samples from the pilot had fiinished by 4:48 pm
#2 samples from the pilot were still running at 8:26 am: 38958375_41,  38958375_47 
#38958375_47  was still running at 10:48 am (21h)
#note: assemblies generally all worked, including co-assemblies within subject, within site
#assemblies generally finished within 24h
#suggested config: --mem=128g --cpus-per-task 8 --time 36:00:00
#################################################### 
cd /data/$USER/metagenomes_mags

################2a. create spades.sh
mkdir /data/$USER/metagenomes_mags/data/01_assembly_single
cat << EOF > spades.sh
	#!/usr/bin/bash
	cd /data/$USER/metagenomes_mags/data/01_assembly_single
	READ1=\$1 
	READ2=\$2
	OUTPUT=\$3
	module load spades
	spades.py --meta -1 \${READ1} -2 \${READ2} -o \${OUTPUT} -t 36 -m 500
EOF

################2b. create megahit.sh
cat << EOF > megahit.sh
	#!/usr/bin/bash
	cd /data/$USER/metagenomes_mags/data/01_assembly_single
	READ1=\$1 
	READ2=\$2
	OUTPUT=\$3
	module load megahit
	megahit -1 \${READ1} -2 \${READ2} -o \${OUTPUT} --memory 128 -t 16 
EOF


	
################2c. launch spades.sh and/or megahit.sh as a swarm 
#create a list of forward reads, reverse reads, and spades_out files
cd /data/$USER/metagenomes_mags
ls -d /data/$USER/metagenomes_mags/data/00_reads/single/*_1.fastq > READ1.list
ls -d /data/$USER/metagenomes_mags/data/00_reads/single/*_2.fastq > READ2.list

#read the files we just created into two arrays
mapfile -t READ1 < READ1.list
mapfile -t READ2 < READ2.list

#let's remove the _1.fastq and replace it with _megahit_out and create an array for the output directories
cut -d/ -f 6 READ1.list > samples
sed 's/_1.fastq/_megahit_out/' samples > MEGAHIT_OUT.list
mapfile -t MEGAHIT_OUT < MEGAHIT_OUT.list
rm samples


#create a megahit swarm file
total=${#READ1[*]}
for (( i=0; i<=$(( $total -1 )); i++ ))
do
echo "bash megahit.sh ${READ1[i]}  ${READ2[i]} ${MEGAHIT_OUT[i]} --presets meta-large --min-contig-len 500 --exclude=cn0863 --memory 128  -t "\$SLURM_CPUS_PER_TASK" --tmp-dir "/lscratch/\$SLURM_JOB_ID""
done > megahit.swarm


#create a spades swarm file
#loop over all the files 
total=${#READ1[*]}
for (( i=0; i<=$(( $total -1 )); i++ ))
do
echo "spades.py --meta -1 ${r1_input[i]} -2 ${r2_input[i]} -o ${spades_out[i]}  -t 8 -m 128"
done > spades.swarm



######################################
#######2d. #launch the swarm job 
#check to make sure the file "spades.swarm" looks good and launch it if so
swarm -f spades.swarm --job-name spades -t 36 -g 800 --module spades --time 36:00:00 --logdir spades_swarm_out -b 3 --partition largemem
swarm -f megahit.swarm --job-name megahit -t 16 -g 128 --module megahit --time 2-0  --logdir megahit_swarm_out_gres_128  --gres=lscratch:400 

######## step 2e.  let's look at the log files
#look at the megahit logs
cd /data/$USER/metagenomes_mags/data/01_assembly_single
ls -d /data/$USER/metagenomes_mags/data/01_assembly_single/**/log > megahit.slogs
mapfile -t logs < megahit.slogs
total=${#logs[*]}
for (( i=0; i<=$(( $total -1 )); i++ ))
do
	awk '/ALL/&&/DONE/' IGNORECASE=1 ${logs[i]}
done > metahit_finished.txt
cat metahit_finished.txt

#look at the spades logs
cd /data/$USER/metagenomes_mags/data/01_assembly_single
ls -d /data/$USER/metagenomes_mags/data/01_assembly_single/**/spades.log > spades.logs
mapfile -t logs < spades.logs
total=${#logs[*]}
for (( i=0; i<=$(( $total -1 )); i++ ))
do
	echo "output from: ${logs[i]}"
	awk '/SPAdes/&&/pipeline/&&/finished/' IGNORECASE=1 ${logs[i]}
done > myspadeslogs_finished.txt
cat myspadeslogs_finished.txt


########step 2f. look at the spades warning files 
#for single assemblies
cd /data/$USER/metagenomes_mags
cat /data/$USER/metagenomes_mags/data/01_assembly_single/**/warnings.log > myspadeswarnings.txt

########2g. get a count of the contigs within each assembly - megahit
#for single assemblies
#Assembled contigs are in /data/$USER/metagenomes_mags/data/01_assembly_single/**/contigs.fasta
ls /data/$USER/metagenomes_mags/data/01_assembly_single/*_megahit_out/final.contigs.fa | wc -l #54 files

#for each of the contigs.fasta files count the number of > within each file
ls -d /data/$USER/metagenomes_mags/data/01_assembly_single/*_megahit_out/final.contigs.fa > mycontigs_megahit.list

#for each of the contigs.fasta files count the number of > within each file
mapfile -t contigs < /data/$USER/metagenomes_mags/data/mycontigs_megahit.list
total=${#contigs[*]}
for ((i=0; i<=$(($total -1)); i++))
do
	egrep -c ">" ${contigs[*]}
done > number_of_contigs_megahit.txt


#################################################### 
#step 3: binning
#note: binning requires a bit more memory for some samples, including the co-assemblies
#8/54 pilot samples failed with 96g ram on 8 cpus and the process seems to be quick, within a couple hours
#the co-assemblies generally worked on 300 gb ram and 16 threads
#suggested config: --mem=356 --cpus-per-task 16 --time 8:00:00
#################################################### 
sinteractive --cpus-per-task=8 --time 36:00:00  --mem=128g --gres=lscratch:400 

#download the metawrap container from sara's paper: https://pubmed.ncbi.nlm.nih.gov/33864056/
module load singularity
#Loading singularity  3.8.5-1  
	#specify the filepath where singularity containers for this project will be housed
	container=/data/$USER/metagenomes_mags/container
	mkdir $container
	cd $container
	#note: you only need to run this once; 
	singularity pull shub://sskashaf/MAG_wf_containers:metawrap


#####3a. define inputs and outputs
#create assembly list which we will use as input along with READ1.list and READ2.list from above
cd /data/$USER/metagenomes_mags/data
	cat /data/$USER/metagenomes_mags/data/mycontigs_megahit.list|wc -l #54

#make an output directory
bins_path=/data/$USER/metagenomes_mags/data/02_binning
	mkdir $bins_path

#####3b. specify inputs and outputs -- single assemblies
cd /data/$USER/metagenomes_mags/data/02_binning
cat << EOF > binning.sh
	#!/usr/bin/bash
	ASSEMBLY=\$1
	OUTPUT=\$2
	READ1=\$3 
	READ2=\$4
	module load singularity
	source /usr/local/current/singularity/app_conf/sing_binds
	singularity run /data/$USER/metagenomes_mags/container/MAG_wf_containers_metawrap.sif \
	metawrap binning  --metabat2 --maxbin2 --concoct -l 5000 -t 32 -m 128 -a \${ASSEMBLY} -o \${OUTPUT} \${READ1} \${READ2}
EOF


#####3c. generate and run the swarm file
#get an array of all the forwad reads in the present working directory
mapfile -t read1 < READ1.list
mapfile -t read2 < READ2.list
mapfile -t assembly < mycontigs_megahit.list  
cut -d/ -f 6 READ1.list > samples
sed 's/_1.fastq/_megahit_binning_out/' samples > BINNING_OUT.list
mapfile -t binning_out < BINNING_OUT.list
rm samples


#write the swarm file
total=${#read1[*]}
for (( i=0; i<=$(( $total -1 )); i++ ))
do
echo "bash binning.sh ${assembly[i]} ${binning_out[i]} ${read1[i]} ${read2[i]}"
done > binning.swarm

#launch the swarm file
swarm -f binning.swarm --job-name binning -t 36 -g 75 --time 30:00:00 --logdir binning_out


################################
#####3d. check that the binning process completed  by examining the output files
cd /data/$USER/metagenomes_mags/data/02_binning/binning_out
ls -d /data/$USER/metagenomes_mags/data/02_binning/binning_out/*o > binning.slurm.out.list
mapfile -t logs < binning.slurm.out.list
total=${#logs[*]}
for (( i=0; i<=$(( $total -1 )); i++ ))
do
	echo "output from: ${logs[i]}"
	awk '/BINNING/&&/PIPELINE/&&/SUCCESSFULLY/&&/FINISHED/' IGNORECASE=1 ${logs[i]}
done > MYBINNING.LOGS
cat  MYBINNING.LOGS


#################################################### 
#step 4: bin refinement
#suggested config: --mem=128g --cpus-per-task 32 --time 12:00:00
#################################################### 
sinteractive --cpus-per-task=8 --time 36:00:00  --mem=128g --gres=lscratch:400 
#make an output directory
refined_bin_path=/data/$USER/metagenomes_mags/data/03_refinement
mkdir $refined_bin_path

cd /data/$USER/metagenomes_mags/data/03_refinement

#make the refinement script
cat << EOF > refinement.sh
	#!/usr/bin/bash
	REFINE_OUT=\$1
	metabat2_bins=\$2
	maxbin2_bins=\$3
	concoct_bins=\$4 
	module load singularity
	source /usr/local/current/singularity/app_conf/sing_binds
	mkdir \${REFINE_OUT} && cd \$_ 
	singularity run /data/$USER/metagenomes_mags/container/MAG_wf_containers_metawrap.sif \
		metawrap bin_refinement -t 32 -o ./ \
		-A  \${metabat2_bins} -B  \${maxbin2_bins} -C  \${concoct_bins}  -c 50 -x 10
EOF


#####4a. launch refinement.sh as swarm
#get the input files into arrays
#we just need the bins, the output from metawrap binning
ls -d /data/$USER/metagenomes_mags/data/02_binning/**/metabat2_bins > metabat2_bins.list 	#use vim to remove the failures dir
ls -d /data/$USER/metagenomes_mags/data/02_binning/**/maxbin2_bins > maxbin2_bins.list  	#use vim to remove the failures dir
ls -d /data/$USER/metagenomes_mags/data/02_binning/**/concoct_bins > concoct_bins.list  	#use vim to remove the failures dir
cut -d/ -f 6 concoct_bins.list > samples
sed 's/_megahit_binning_out/_refinement_out/' samples > REFINEMENT_OUT.list
rm /data/$USER/metagenomes_mags/data/samples


#note that these 3 lists should be equal in length and they should have the same samples in the same order across binners
cat metabat2_bins.list |wc -l #194
cat maxbin2_bins.list |wc -l #194
cat concoct_bins.list |wc -l #194
cat REFINEMENT_OUT.list|wc -l #194
 
#####4b. generate swarm file for refinement
mapfile -t refine_out < REFINEMENT_OUT.list
mapfile -t metabat2_bins < metabat2_bins.list
mapfile -t maxbin2_bins < maxbin2_bins.list 
mapfile -t concoct_bins < concoct_bins.list


#launch swarm
total=${#refine_out[*]}
for (( i=0; i<=$(( $total -1 )); i++ ))
do
echo "bash refinement.sh ${refine_out[i]} ${metabat2_bins[i]} ${maxbin2_bins[i]} ${concoct_bins[i]}"
done > refinement.swarm
swarm -f refinement.swarm --job-name refines -t 32 -g 128 --time 24:00:00 --logdir refinement_out




#### 4c. check how many refined bins we have
cat /data/$USER/metagenomes_mags/data/03_refinement/**/metawrap_50_10_bins.stats| awk '$2>50 && $3<10' | wc -l #213
cat /data/$USER/metagenomes_mags/data/03_refinement/**/metawrap_50_10_bins.stats| awk '$2>75 && $3<5' | wc -l #187

### 4d. get a list of bins
ls -d /data/$USER/metagenomes_mags/data/03_refinement/**/metawrap_50_10_bins/*fa > refined_bins.list #179


#################################################### 
#step 5: QC MAGs (GUNC, infernal, trnascan, checkm2), rename them, and dereplicate
#warning: a manual step was done. this requires caution.
#I cheated here and added some columns to the file bins_refined.lists manually in Excel
#importantly, i separated the path so that the desired BIN name is identified in column 5 of the files: 1) bins_refined.list.txt and 2) bins_refined.coassembly.list.txt
#this will allow us to copy the bin from the original path to files identified by the array ${NewNAME[*]}
#suggested config: sinteractive --cpus-per-task=8 --time 36:00:00  --mem=128g --gres=lscratch:400 
#################################################### 
#####step 5a: create a drep conda environment and add the checkM data path
#note: you only have to do this once unless you move the checkm database
sinteractive --cpus-per-task=8 --time 36:00:00  --mem=128g --gres=lscratch:400 
conda create -n drep2
conda activate drep2
conda config --add channels bioconda; conda install drep
export CHECKM_DATA_PATH=/data/$USER/metagenomes_mags/data/databases/MY_CHECKM_FOLDER



#####5b. single bins-rename bins and copy them into input directory, one central place 
derep_bin_path=/data/$USER/metagenomes_mags/data/04_dereplication
mkdir $derep_bin_path
mkdir $derep_bin_path/input

cd /data/$USER/metagenomes_mags/data

#because i cheated and added manual rows we need to remove special characters that excel adds to file, specifically ^M
sed -e "s/\r//g" bins_refined.list > bins_refined.list.txt
MYBINPATH=($(awk '{print $1}' bins_refined.list.txt)) #this should be a column containing the path to the bin
NewNAME=($(awk '{print $5}' bins_refined.list.txt)) #this should be the path with / replaced by _ so that we have the sample name with the bin name

total=${#MYBINPATH[*]}
for (( i=0; i<=$(( $total -1 )); i++ ))
do
	echo "cp ${MYBINPATH[i]} /data/$USER/metagenomes_mags/data/04_dereplication_all/input/${NewNAME[i]}"
done  > renamebins.sh
bash renamebins.sh



#####5c. Run GUNC on prokaryotic MAgs. 
#we copy the names of these mags into a file called gunc_passes.txt
cat << EOF > gunc.sh
#!/usr/bin/bash
source myconda
conda activate gunc
gunc run --input_dir /data/$USER/metagenomes_mags/data/04_dereplication/input/ -r /data/Segrelab/Diana/databases/gunc/gunc_db_progenomes2.1.dmnd --out_dir ../gunc_out;
EOF
sbatch --mem=128g --cpus-per-task 30 --time 36:00:00 gunc.sh 


#####5d. Run checkM2 on prokaryotic MAGs.
#https://github.com/chklovski/CheckM2
cat << EOF > checkm2.sh
#!/usr/bin/bash
source myconda
conda activate checkm2
/data/$USER/metagenomes_mags/checkm2/bin/checkm2 predict --threads 30 --input /data/$USER/metagenomes_mags/data/04_dereplication/input --output-directory ../checkm2_out -x fa
EOF
sbatch --cpus-per-task 50 --mem=500g --time 8:00:00 --partition largemem checkm2.sh --job-name checkm


#####5e. Run GUNC.rmd to generate a file of GUNC and checkm2 passes
mkdir /data/$USER/metagenomes_mags/data/04_dereplication/qc_pass

#####5f. Collect the  MAGs that pass GUNC and checkM2
passes.list = QC_pass_guncAndCheckM2.tsv 
mapfile -t passed < passes.list
total=${#passed[*]}
for (( i=0; i<=$(( $total -1 )); i++ ))
do
echo "cp input/${passed[i]} /data/$USER/metagenomes_mags/data/04_dereplication/qc_pass"
done > getqcpasses.sh #1438 lines
bash getqcpasses.sh
 

#####5g. dereplicate prokaryotic MAgs that pass both GUNC and checkm2
#there are 1438 that do
#cd /data/$USER/metagenomes_mags/data/04_dereplication
cat << EOF > derep.sh
	#!/usr/bin/bash
	source myconda
	conda activate drep
	dRep dereplicate -p 32 /data/$USER/metagenomes_mags/data/04_dereplication/derep_out/ -g /data/$USER/metagenomes_mags/data/04_dereplication/qc_pass/*fa -pa 0.90 -sa 0.95 -nc 0.30 -cm larger -comp 50 -con 10
EOF
sbatch --mem=128g --cpus-per-task 32 --time 3-0 derep.sh


#####5h. Run Infernal and tRNA-SCAN to get rRNA info

#cd /data/$USER/metagenomes_mags/data/databases
#wget ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.cm.gz
#gunzip Rfam.cm.gz
#wget ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.clanin
#cmpress Rfam.cm

cd /data/$USER/metagenomes_mags/data/04_dereplication/input
mybins=($(ls  *fa)) 
tnra_names=($(echo "${mybins[*]}" | sed 's/.fa/_trnascan.txt/g'))
infernal_names=($(echo "${mybins[*]}" | sed 's/.fa/_infernal_out.txt/g'))
subtotal=${#mybins[*]}


#### Run infernal
source myconda
conda activate infernal
for (( j=0; j<=$(( $subtotal -1 )); j++));
do echo "cmsearch -Z 1000 --hmmonly --cut_ga --noali --tblout=${infernal_names[j]} /data/$USER/metagenomes_mags/data/databases/Rfam.cm input/${mybins[j]};"
done > infernal.sh; 
nano infernal.sh #add shebang and activate conda env

#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6768409/
### Run trna-scan
conda activate trnascan
for (( j=0; j<=$(( $subtotal -1 )); j++));
do 
echo "tRNAscan-SE  input/${mybins[j]} -B -Q -o ${tnra_names[j]}";
done > transcan.sh;
nano transcan.sh #add shebang and activate conda envs
sbatch --mem=128g --cpus-per-task 8 --timme 36:00:00 transcan.sh



#################################################### 
#step 6: taxonomic classification of prokaryotes: https://ecogenomics.github.io/GTDBTk/installing/bioconda.html 
#suggested config: for 504 MAGS: sbatch --mem=128g --cpus-per-task 8 --time 36:00:00 
#note the output file genomeInformation.csv will be needed for the script: mapping_bwa2abundance.rmd
#################################################### 

######6a. run gtdbk on the full collection of mags, irresspective of whether they passed gunc or checkm2 
#get version 2 of gtdbk which is required for database 207
source myconda
#note that gatdk version  must be compatible with the database verrsion or 
#conda create -n gtdbtk-2.0.0 -c conda-forge -c bioconda gtdbtk=2.0.0
conda activate gtdbtk-2.0.0 

######6b. download database, release 207
#see: https://data.ace.uq.edu.au/public/gtdb/data/releases/
GTDBTK_DB_PATH=/data/$USER/conda/envs/gtdbtk-2.0.0/share/gtdbtk-2.0.0/db/
cd $GTDBTK_DB_PATH

###### GTDBTK_DB_PATH is defined in build.sh, store the db where it's expected
wget https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/auxillary_files/gtdbtk_data.tar.gz
tar xvzf ${GTDBTK_DB_PATH}/gtdbtk_data.tar.gz -C ${GTDBTK_DB_PATH} --strip 1

######6c. set up
#the dereplicated bins live here: /data/$USER/metagenomes_mags/data/04_dereplication_all/derep_out/dereplicated_genomes

######6d.  run gatkdb on ALL the assemblies that pass QC
mkdir /data/$USER/metagenomes_mags/data/05_gatkdb

cat << EOF > gtdbk.sh
#!/usr/bin/bash
source myconda
conda activate gtdbtk-2.0.0 
GTDBTK_DB_PATH=/data/$USER/conda/envs/gtdbtk-2.0.0/share/gtdbtk-2.0.0/db/
input_genomes=/data/$USER/metagenomes_mags/data/04_dereplication/qc_pass
output_directory=/data/$USER/metagenomes_mags/data/05_gatkdb
gtdbtk classify_wf --cpus 16 --genome_dir  \$input_genomes --out_dir  \$output_directory -x fa
EOF
sbatch --mem=128g --cpus-per-task 16 --time 36:00:00  gtdbk.sh 




#################################################### 
#step 7: taxonomic classification of eukaryotes: eukcc on individual directories
#https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02155-4#Sec13
#suggested config:  sbatch --mem=128g --cpus-per-task 8 --time 6:00:00 "$f" 
#################################################### 
sinteractive --mem=128g --cpus-per-task 8 --time 36:00:00


######7a. download eukcc and specify the path to the database
#EukCC version 2.1.0
cd /data/$USER/metagenomes_mags/data/databases
mkdir eukccdb
wget http://ftp.ebi.ac.uk/pub/databases/metagenomics/eukcc/eukcc2_db_ver_1.1.tar.gz
tar -xzvf eukcc2_db_ver_1.1.tar.gz
export EUKCC2_DB=/data/$USER/metagenomes_mags/data/databases/eukcc2_db_ver_1.1  
#We assume that you did set you $EUKCC2_DB to the correct location. If not please use the --db flag to pass the database to EukCC.

######7b. prepare to write a bunch of scripts
cd /data/proctordm/mags_abx_batch/04_derep/06_eukcc


mapfile -t bins <bins.list
mapfile -t output < output.list


######7c. generate a script to run eukcc on each sample
cat << EOF > eukcc.sh
#!/usr/bin/bash
eukcc_out=\$1
FOLDER=\$2
source myconda
conda activate eukcc
export EUKCC2_DB=/data/$USER/metagenomes_mags/data/databases/eukcc2_db_ver_1.1  
mkdir \${eukcc_out} 
cd \${eukcc_out}
eukcc folder \${FOLDER} --out /data/$USER/metagenomes_mags/data/06_eukcc/\${eukcc_out} --threads 16
EOF



######7e. download eukcc and specify the path to the database
#EukCC version 2.1.0
cd /data/$USER/metagenomes_mags/data/databases
mkdir eukccdb
wget http://ftp.ebi.ac.uk/pub/databases/metagenomics/eukcc/eukcc2_db_ver_1.1.tar.gz
tar -xzvf eukcc2_db_ver_1.1.tar.gz
export EUKCC2_DB=/data/$USER/metagenomes_mags/data/databases/eukcc2_db_ver_1.1  
#We assume that you did set you $EUKCC2_DB to the correct location. If not please use the --db flag to pass the database to EukCC.


##6c. generate a script to run eukcc on each CONCOCT BIB
mkdir /data/$USER/metagenomes_mags/data/06_eukcc

cat << EOF > eukcc.sh
#!/usr/bin/bash
eukcc_out=\$1
FOLDER=\$2
source myconda
conda activate eukcc
export EUKCC2_DB=/data/$USER/metagenomes_mags/data/databases/eukcc2_db_ver_1.1  
mkdir \${eukcc_out} && cd \$_ 
eukcc folder \${FOLDER} --out /data/$USER/metagenomes_mags/data/06_eukcc/euk_pass/\${eukcc_out} --threads 16
EOF

#lit all the CONCOCT BINS
mapfile -t contigs <   concoctbins.txt
#define outputs
output=($(echo "${contigs[*]}" | sed 's/.fa/_eukcc_out/g'))


#create a eukcc swarm job, running it on all the CONCOCT BINS
total=${#contigs[*]}
for ((i=0; i<=$(($total -1)); i++))
do
	echo "bash eukcc.sh  ${output[i]}  ${contigs[i]}"
done > eukcc.swarm
swarm -f eukcc.swarm  --job-name eukcc -t 16 -g 128 --time 24:00:00 -b 2 --logdir eukcc_out


#####6e. add the directory name to each of the important output files, eukcc.csv
target=/data/$USER/metagenomes_mags/data/06_eukcc/eukcc_csv
mkdir $target
cd $target

#get all the eukaryotic csv files
ls -d /data/$USER/metagenomes_mags/data/06_eukcc/**/*eukcc.csv > euk.csv.files

#get sammple names so we can rename csv files with their sample name and add a line with the file path to the directory
cat euk.csv.files | cut -d/ -f7  > samples
mapfile -t samples < samples
mapfile -t FILE <  euk.csv.files
NEWFILE=($(echo "${samples[*]}" | sed 's/_eukcc_out/.rev.eukcc.csv/g'))


#now let's loop over the files and 1) rename files to NEWFILE, 2) add filepath PATH and 3) copy the revised file NEWFILE to $target
total=${#samples[*]}
for ((i=0; i<=$(($total -1)); i++))
do
	sed 's,^,'${FILE[i]}' ,' ${FILE[i]} > ${NEWFILE[i]}
done #now e have the csv files in one place

#concatenate the .csv files dropping the first line
tail -n +2 -q /data/$USER/metagenomes_mags/data/06_eukcc/eukcc_csv/*csv >> euk.bins.csv #inspect this

#6f. get the eukaryotic bins in one place 
bin_destiny=/data/$USER/metagenomes_mags/data/06_eukcc/input
mkdir $bin_destiny
cd $bin_destiny

#input euk bins
ls -d /data/$USER/metagenomes_mags/data/02_binning/**/concoct_bins/*fa > concoctbins.txt


mapfile -t BINS <   concoctbins.txt
total=${#BINS[*]}
for ((i=0; i<=$(($total -1)); i++))
do
	path=${BINS[i]}
	foo=${path//\//_}
	NEWNAME=($(echo "$foo" | sed 's/_data_proctordm_metagenomes_mags_data_02_binning_//g'))
	cp ${BINS[i]} ./$NEWNAME
done #this contains all the concoct bins



#6g. In excel, open the file: eukcc_csv/euk.bins.csv
#copy the path into column B
#separate column B into components
#use concatenate to create a column in the euk.bins.csv file with filenames in the format of $NEWNAME 'Met4805_bins_concoct_bins_bin.3.fa'
#subset on bins >50% completeness and < 10% contamination. We shouold have 107 bins iin a new folder
#save the file as a new text file: euk.bins_pass.txt 
#change $17 to the column number where $NEWNAME is located
awk '{ print $7 }' /data/$USER/metagenomes_mags/data/06_eukcc/eukcc_csv/euk.bins.csv > bins2keep
sed 's/\r$//' bins2keep > bins2keep2 #remove special characters introduced by excel

ALLBINS=data/$USER/metagenomes_mags/data/06_eukcc/input
SUBSET=data/$USER/metagenomes_mags/data/06_eukcc/qc_pass
mkdir $SUBSET
cd $ALLBINS

mapfile -t BINS <   ../bins2keep2
total=${#BINS[*]}
for ((i=0; i<=$(($total -1)); i++))
do
	cp "${BINS[i]}" $SUBSET
done

cd $SUBSET


#note all 107 eukaryotic bins live here
#/data/$USER/metagenomes_mags/data/06_eukcc/eukcc_pass

#################################################### 
#step 8: taxonomic classification of eukaryotes: mash  
#note; you must use scratch space since we're launching an r job with bash
#################################################### 
sinteractive --mem=96g --cpus-per-task 8 --time 36:00:00 --gres=lscratch:400
mkdir /data/$USER/metagenomes_mags/data/08_mash_euk
cd /data/$USER/metagenomes_mags/data/08_mash_euk


##### step 8a: run mash against refseq or genbank
module load mash
#mash against genbank
 
mysketch=/data/$USER/metagenomes_mags/data/databases/fungi.genbank.msh
mash dist -p 8 $mysketch /data/$USER/metagenomes_mags/data/06_eukcc/eukcc_pass/*fa > fungi.genbank_mash_out.txt


##### step 8b: make a bash script to run the R script on HPC
#note: i had to run this on the command line, not with sbatch or bash 
cat << EOF > launch_getBestMashHit.sh
#!/bin/bash
module load R
R --vanilla < ./getBestMashHit.r > ./MashR.out
EOF

##### step 8c: make an r script to get best mash hit for each MAG (we just minimize p. User wants to look at hash matches / mash distance as well)
cat << EOF >  getBestMashHit.r 
library("data.table")
OMP_NUM_THREADS=8


dat=data.table::fread("fungi.genbank_mash_out.txt")
colnames(dat) = c("Reference", "MAG", "distance", "p.value", "hash")

mymags=unique(dat$MAG)
holder <- vector("list", length(mymags))
names(holder) = mymags
  for(i in 1:length(mymags)) {
 	df=subset(dat, MAG==mymags[i])
 	df = df[which.min(df$p.value),]
 	holder[[i]] = df
 	}

bestmashhits=do.call("rbind", holder)
bestmashhits$Line = paste0("Line", 1:nrow(bestmashhits))
write.table(bestmashhits, file="bestmashhits.fungi.genbank.msh_manual.txt")
EOF


##### step 8d: run the r script
bash launch_getBestMashHit.sh


##### step 8e: fix the output file by removing the header row
tail -n +2 bestmashhits.fungi.genbank.msh_manual.txt> bestmashhits.fungi.genbank.msh_noheader.txt
sed 's/\"//g' bestmashhits.fungi.genbank.msh_noheader.txt  > bestmashhits.fungi.genbank.msh_noheader_noquotes.txt

#################################################### 
#step 9: taxonomic classification of eukaryotes: mummer  
#note the output file fungal_taxonomy.csv will be needed for the script: mapping_bwa2abundance.rmd
#sara's recommendation is to use > 30% aligned fraction;> 95% ANI to identify the taxon
#note: we will lose taxa in our final table that don't appear in the file: assembly_summary_genbank_fungi.txt
#this meeans any taxa from smgc and other misc genomes need to be added back
#be sure to keep track of how many mags you start with and how many you get taxonomic classifications for
# we can prob fix this by using plyr::full_join in R, rather than merge
#################################################### 
cd /data/$USER/metagenomes_mags/data/08_mash_euk

module load mummer
#[+] Loading mummer  4.0.0beta2  on cn0920 

##### step 9a: write a bunch of mummer scripts we can launch

cat << EOF > mummer.sh
#!/usr/bin/bash
REFERENCE=\$1 #positional parameters; this wil allow you to enter the arguments from 1:4
MAG=\$2
OUT=\$3
module load mummer
mkdir \${OUT}_mummer_output && cd \$_ #the $_ specifies cd into the previous argument
dnadiff \${REFERENCE} \${MAG} -p \${OUT}
file=\${OUT}.rev
cat *report | grep -i "^AlignedBases" > \$file
cat *report | grep -i "^AvgIdentity" | awk 'NR%2!=0'>>  \$file 
cp \$file /data/$USER/metagenomes_mags/data/08_mash_euk
EOF
#the backslashes here will tell it to create the file 
#create a swarm for each element of REFERENCE, MAG, OUT, target

#set ourselves up with a bunch of arrays by reading in the output of the R script which got the lowest p value per MAG
input=bestmashhits.fungi.genbank.msh_noheader_noquotes.txt


#load the reference genomes fosr p < 0.05 into an array 
REFERENCE=($(awk '{print $2}' $input))
#load the MAGS for p < 0.05 into an array s
MAG=($(awk '{print $3}' $input))
#output
OUT=($(awk '{print $7}' $input)) 

total=${#REFERENCE[*]}
for ((i=0; i<=$(($total -1)); i++))
do
	echo "bash mummer.sh ${REFERENCE[i]} ${MAG[i]} ${OUT[i]}"
done > mummer.swarm

swarm -f  mummer.swarm  --job-name mummer -t 16 -g 98 --time 2:00:00 -b 2 --logdir mummer_swarm_out


##### step 9c: aggregate the mummer reports in R
cd /data/$USER/metagenomes_mags/data/08_mash_euk
module load R
R
library(data.table)
path <- "/data/$USER/metagenomes_mags/data/08_mash_euk" 

# Forward and reverse fastq filenames have format: cat.tsv, the output of the blast searches on the 4 genomes
mummerOut <- sort(list.files(path, pattern=".rev", full.names = TRUE))


#read in the files and store them in a list
TSV.list <- vector("list", length(mummerOut))

for(i in 1:length(mummerOut)) {
      myDat = data.frame(read.table(mummerOut[[i]], header = FALSE))
      myDat = data.frame(t(myDat))
      colnames(myDat) <- myDat[1,]
      myDat <- myDat[-1,]
      myDat$mummerLink = basename(mummerOut[[i]])
      TSV.list[[i]] = myDat
}

mummerDF = do.call("rbind", TSV.list)
mummerDF$mummerLink = stringr::str_remove_all(mummerDF$mummerLink, ".rev")
mashTable = read.table("/data/$USER/metagenomes_mags/data/08_mash_euk/bestmashhits.fungi.genbank.msh_noheader_noquotes.txt")
mashTable$Reference = basename(mashTable$V2)
mashTable$MAG = basename(mashTable$V3)
colnames(mashTable) = c("Junk","ReferencePath", "MAGPath", "Distance", "p.value", "Hash","mummerLink", "Reference", "MAG")

myFungi = merge(mashTable, mummerDF, by="mummerLink")
myFungi$Reference2 = stringr::str_remove_all(myFungi$Reference, "_genomic.fna")
ncbi.tax.table= data.table::fread("/data/$USER/metagenomes_mags/data/databases/assembly_summary_genbank_fungi.txt")

myFungi3 = merge(myFungi, ncbi.tax.table, by="Reference2")
write.csv(myFungi3, file="fungal_taxonomy_v2.0.csv")
q()

   
#################################################### 
#step 9: aggregate the MAGs AND fix the headers
#if you don't fix thee headers you will get bwa mapping results to each contig within bins in the subsequent steps
#by fixing the headers you get mapping results to bins in later steps
#################################################### 
#######step 10a: get the mags in one place, combine eukaryotic mags with prokaryotic ones
mkdir /data/$USER/metagenomes_mags/data/09_mapping
cp /data/$USER/metagenomes_mags/data/04_dereplication/derep_out/* 09_mapping/
cp /data/$USER/metagenomes_mags/data/07_euk_derep/* 09_mapping/
ls *fa |wc -l

cd /data/$USER/metagenomes_mags/data/09_mapping
mkdir header


#######step 10B; fix the bin names so that reads map to bins instead of contigs
#fix eukaryotic bin names
sourcemyconda 
conda activate py2
inp=`ls | grep fa""`
for BIN in $inp;
do python /home/$USER/bin/MGS-gut/scripts/rename_multifasta_prefix.py -f $BIN -p $BIN > header/$BIN;
done
cd header
cat *fa > ../eukBacBins.derep.fa
mkdir input
mv Met*fa input




################################################### 
#step 12: map reads against the dereplicated MAGS
#we use bwa to make reads to the mag collection
#This chunk relies heavily on Ryan Blaustein's code; I just turned it into a script that was executable by swarm
# we reply on two scripts from Alex Alemeida's github map2ref.sh and parse_bwa-depth.py
#################################################### 
cd /data/$USER/metagenomes_mags/data/09_mapping



########12a. index the reference
#index the bins with bwa
module load bwa
#[+] Loading bwa 0.7.17 
bwa index eukBacBins.derep.fa




#######step 12b create a copy of Alex Almeida's  script modified to run on biowulf
#need to use a conda env running python 2.7 with numpy installed

cat << EOF > map.sh
#!/bin/bash
module load bwa
module load samtools
source myconda
conda activate py2

if [ $# -eq 0 ]; then
    echo "Map reads against a custom reference database"
    echo ""
    echo "Notes:" 
    echo "- run 'bwa index' on reference FASTA beforehand"
    echo "- headers in FASTA file must include genome name (use rename_multifasta_prefix.py script)"
    echo "" 
    echo "usage: script.sh input_1.fastq(gz) input_2.fastq(gz) ref.fasta outprefix"
    echo ""
    exit 1
fi

# variables
ref=${3}
refix=$(basename ${ref%%.fa*})
reads=${1}
reads2=${2}
outprefix=${4}
readname=$(basename ${outprefix})

# initial mapping
bwa mem -t 8 ${ref} ${reads} ${reads2} | samtools view -@ 7 -uS - -o ${outprefix}_${refix}_unsorted.bam

# sort bam file
samtools sort -@ 7 ${outprefix}_${refix}_unsorted.bam -o ${outprefix}_${refix}_sorted.bam

# extract unique counts
samtools view -@ 7 -q 1 -f 2 -u ${outprefix}_${refix}_sorted.bam -o ${outprefix}_${refix}_unique_sorted.bam
samtools index -@ 7 ${outprefix}_${refix}_unique_sorted.bam
samtools idxstats ${outprefix}_${refix}_unique_sorted.bam > ${outprefix}_${refix}_unique_depth.tab; samtools depth ${outprefix}_${refix}_unique_sorted.bam > ${outprefix}_${refix}_unique_depth-pos.tab; 
python parse_bwa-depth.py ${outprefix}_${refix}_unique_depth.tab ${outprefix}_${refix}_unique_depth-pos.tab > ${outprefix}_${refix}_unique.tab; rm -rf ${outprefix}_${refix}_unique_sorted.ba* ${outprefix}_${refix}_unique_depth*

# extract total counts
samtools index -@ 7 ${outprefix}_${refix}_sorted.bam
samtools idxstats ${outprefix}_${refix}_sorted.bam > ${outprefix}_${refix}_depth.tab; samtools depth ${outprefix}_${refix}_sorted.bam > ${outprefix}_${refix}_depth-pos.tab; 
python parse_bwa-depth.py ${outprefix}_${refix}_depth.tab ${outprefix}_${refix}_depth-pos.tab > ${outprefix}_${refix}_total.tab; rm -rf ${outprefix}_${refix}_unsorted.bam ${outprefix}_${refix}_sorted.bam.bai ${outprefix}_${refix}_depth*
EOF


###################12c. make a swarm file to run read mapping for each sample
ls -d /data/$USER/metagenomes_mags/data/00_reads/single/*_1.fastq > READ1.list
ls -d /data/$USER/metagenomes_mags/data/00_reads/single/*_2.fastq > READ2.list

mapfile -t READ1 < READ1.list
mapfile -t READ2 < READ2.list


#let's set ourselves up with script names for each sample
cut -d/ -f6 READ1.list | sed 's/_1.fastq/_mapping/g' > samples
mapfile -t samples < samples

total=${#READ1[*]}
for ((i=0; i<=$(($total -1)); i++))
do
	echo "bash map.sh ${READ1[i]} ${READ2[i]} eukBacBins.derep.fa  ${samples[i]}"
done > map_reads.swarm
swarm -f map_reads.swarm --job-name mapping -t 20 -g 128 --time 8:00:00 --logdir mapping_swarm_out --gres=lscratch:100 


################################################### 
#step 12:  run srst2 on each sample using the card database
#################################################### 

##6c. 
cat << EOF > srst2.sh
#!/usr/bin/bash
source myconda
conda activate srst2
READ1=\$1
READ2=\$2
FU=\$3
srst2 --input_pe \${READ1} \${READ2} --output \${FU} --log --gene_db CARD_v3.0.8_SRST2.fasta
EOF

mapfile -t read1 < READ1.list
mapfile -t read2 < READ2.list
cat READ1.list | cut -d/ -f6  > samples
mapfile -t samples < samples

output=($(echo "${samples[*]}" | sed 's/_1.fastq/_srst2_out/g'))
total=${#read1[*]}
for ((i=0; i<=$(($total -1)); i++))
do
	echo "bash srst2.sh ${read1[i]} ${read2[i]} ${output[i]}"
done > srst2.swarm

swarm -f srst2.swarm --job-name srst -t 8 -g 18 --time 36:00:00 --logdir srst_swarm_out

################################################### 
#step 13:  run fastani
#################################################### 
ls *fa > query.list
ls *fa > ref.list

#!/usr/bin/bash
source myconda
conda activate fastani
fastANI --ql input/query.list --rl input/ref.list -o output/pairwise_fastani_sbatch.txt


################################################### 
#step 13:  run GotoTree
#################################################### 

#!/usr/bin/bash
source myconda
conda activate gtotree
GToTree -a refseq_accessions.txt \
        -f fasta_files_withPacBio.txt -H Gammaproteobacteria \
        -t -L Species,Strain -m map.tsv -j 32 \
        -o kleb2_gtree_out_otherNH
        
