sbatch  --gres=lscratch:400  --cpus-per-task=2 --mem=128g --time 36:00:00


### see: isolates.csv
### we start with 141 isolates downloaded on Sept. 12/13 20222 from IL, MA, FL
#get the samples that didn't download with SAM numbers
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
#Map reads with bwa
###################################
bwa index reference/Caur.B11891_adqia.v2.Nano.pilon.fasta -p B11891

##make a bwa.sh script that we will run as a swarm
cat << EOF > bwa.sh
	#!/usr/bin/bash
	READ1=\$1 
	READ2=\$2
	OUTPUT=\$3
	module load bwa
	bwa mem -M -t 16 /data/proctordm/cauris_genomes/reference/B11891 \${READ1} \${READ2} > alignments/\${OUTPUT}
EOF
​

#map to SRA
bash bwa.sh /data/proctordm/cauris_genomes/fastq/SRA/SAMN10139260_1.fastq /data/proctordm/cauris_genomes/fastq/SRA/SAMN10139260_2.fastq FL_B12847_urine_2019.sam
bash bwa.sh /data/proctordm/cauris_genomes/fastq/SRA/SAMN16376912_1.fastq /data/proctordm/cauris_genomes/fastq/SRA/SAMN16376912_2.fastq FL_B12847_urine_2021.sam
bash bwa.sh /data/proctordm/cauris_genomes/fastq/SRA/SAMN14925561_1.fastq /data/proctordm/cauris_genomes/fastq/SRA/SAMN14925561_2.fastq FL_B16329_blood_2020.sam
bash bwa.sh /data/proctordm/cauris_genomes/fastq/SRA/SAMN10139552_1.fastq /data/proctordm/cauris_genomes/fastq/SRA/SAMN10139552_2.fastq IL_B11842_urine_2019.sam
bash bwa.sh /data/proctordm/cauris_genomes/fastq/SRA/SAMN10139551_1.fastq /data/proctordm/cauris_genomes/fastq/SRA/SAMN10139551_2.fastq IL_B11843_blood_2019.sam
bash bwa.sh /data/proctordm/cauris_genomes/fastq/SRA/SAMN10139548_1.fastq /data/proctordm/cauris_genomes/fastq/SRA/SAMN10139548_2.fastq IL_B11880_nares_2019.sam
bash bwa.sh /data/proctordm/cauris_genomes/fastq/SRA/SAMN10139547_1.fastq /data/proctordm/cauris_genomes/fastq/SRA/SAMN10139547_2.fastq IL_B11882_nares_2019.sam
bash bwa.sh /data/proctordm/cauris_genomes/fastq/SRA/SAMN10139546_1.fastq /data/proctordm/cauris_genomes/fastq/SRA/SAMN10139546_2.fastq IL_B11883_Rectal_2019.sam
bash bwa.sh /data/proctordm/cauris_genomes/fastq/SRA/SAMN10139545_1.fastq /data/proctordm/cauris_genomes/fastq/SRA/SAMN10139545_2.fastq IL_B11884_Vaginal_2019.sam
bash bwa.sh /data/proctordm/cauris_genomes/fastq/SRA/SAMN10139544_1.fastq /data/proctordm/cauris_genomes/fastq/SRA/SAMN10139544_2.fastq IL_B11885_Axilla_2019.sam
bash bwa.sh /data/proctordm/cauris_genomes/fastq/SRA/SAMN10139543_1.fastq /data/proctordm/cauris_genomes/fastq/SRA/SAMN10139543_2.fastq IL_B11886_Groin_2019.sam
bash bwa.sh /data/proctordm/cauris_genomes/fastq/SRA/SAMN10139542_1.fastq /data/proctordm/cauris_genomes/fastq/SRA/SAMN10139542_2.fastq IL_B12028_blood_2019.sam
bash bwa.sh /data/proctordm/cauris_genomes/fastq/SRA/SAMN10139541_1.fastq /data/proctordm/cauris_genomes/fastq/SRA/SAMN10139541_2.fastq IL_B12029_blood_2019.sam
bash bwa.sh /data/proctordm/cauris_genomes/fastq/SRA/SAMN10139540_1.fastq /data/proctordm/cauris_genomes/fastq/SRA/SAMN10139540_2.fastq IL_B12030_blood_2019.sam
bash bwa.sh /data/proctordm/cauris_genomes/fastq/SRA/SAMN10139539_1.fastq /data/proctordm/cauris_genomes/fastq/SRA/SAMN10139539_2.fastq IL_B12031_blood_2019.sam
bash bwa.sh /data/proctordm/cauris_genomes/fastq/SRA/SAMN10139538_1.fastq /data/proctordm/cauris_genomes/fastq/SRA/SAMN10139538_2.fastq IL_B12032_blood_2019.sam
bash bwa.sh /data/proctordm/cauris_genomes/fastq/SRA/SAMN10139537_1.fastq /data/proctordm/cauris_genomes/fastq/SRA/SAMN10139537_2.fastq IL_B12033_blood_2019.sam
bash bwa.sh /data/proctordm/cauris_genomes/fastq/SRA/SAMN10139536_1.fastq /data/proctordm/cauris_genomes/fastq/SRA/SAMN10139536_2.fastq IL_B12034_blood_2019.sam
bash bwa.sh /data/proctordm/cauris_genomes/fastq/SRA/SAMN10139535_1.fastq /data/proctordm/cauris_genomes/fastq/SRA/SAMN10139535_2.fastq IL_B12035_blood_2019.sam
bash bwa.sh /data/proctordm/cauris_genomes/fastq/SRA/SAMN10139534_1.fastq /data/proctordm/cauris_genomes/fastq/SRA/SAMN10139534_2.fastq IL_B12036_blood_2019.sam
bash bwa.sh /data/proctordm/cauris_genomes/fastq/SRA/SAMN10139527_1.fastq /data/proctordm/cauris_genomes/fastq/SRA/SAMN10139527_2.fastq IL_B12046_AxillaGroin_2019.sam
bash bwa.sh /data/proctordm/cauris_genomes/fastq/SRA/SAMN10139526_1.fastq /data/proctordm/cauris_genomes/fastq/SRA/SAMN10139526_2.fastq IL_B12077_feces_2019.sam
bash bwa.sh /data/proctordm/cauris_genomes/fastq/SRA/SAMN10139513_1.fastq /data/proctordm/cauris_genomes/fastq/SRA/SAMN10139513_2.fastq IL_B12189_nares_2019.sam
bash bwa.sh /data/proctordm/cauris_genomes/fastq/SRA/SAMN10139512_1.fastq /data/proctordm/cauris_genomes/fastq/SRA/SAMN10139512_2.fastq IL_B12190_Axilla_2019.sam
bash bwa.sh /data/proctordm/cauris_genomes/fastq/SRA/SAMN10139511_1.fastq /data/proctordm/cauris_genomes/fastq/SRA/SAMN10139511_2.fastq IL_B12229_kidney_2019.sam
bash bwa.sh /data/proctordm/cauris_genomes/fastq/SRA/SAMN10139501_1.fastq /data/proctordm/cauris_genomes/fastq/SRA/SAMN10139501_2.fastq IL_B12378_AxillaGroin_2019.sam
bash bwa.sh /data/proctordm/cauris_genomes/fastq/SRA/SAMN10139500_1.fastq /data/proctordm/cauris_genomes/fastq/SRA/SAMN10139500_2.fastq IL_B12380_AxillaGroin_2019.sam
bash bwa.sh /data/proctordm/cauris_genomes/fastq/SRA/SAMN10139499_1.fastq /data/proctordm/cauris_genomes/fastq/SRA/SAMN10139499_2.fastq IL_B12388_AxillaGroin_2019.sam
bash bwa.sh /data/proctordm/cauris_genomes/fastq/SRA/SAMN10139498_1.fastq /data/proctordm/cauris_genomes/fastq/SRA/SAMN10139498_2.fastq IL_B12406_urine_2019.sam
bash bwa.sh /data/proctordm/cauris_genomes/fastq/SRA/SAMN10139373_1.fastq /data/proctordm/cauris_genomes/fastq/SRA/SAMN10139373_2.fastq IL_B12510_AxillaGroin_2019.sam
bash bwa.sh /data/proctordm/cauris_genomes/fastq/SRA/SAMN10139351_1.fastq /data/proctordm/cauris_genomes/fastq/SRA/SAMN10139351_2.fastq IL_B12579_Axilla_2019.sam
bash bwa.sh /data/proctordm/cauris_genomes/fastq/SRA/SAMN10139227_1.fastq /data/proctordm/cauris_genomes/fastq/SRA/SAMN10139227_2.fastq IL_B12937_Axilla_2019.sam
bash bwa.sh /data/proctordm/cauris_genomes/fastq/SRA/SAMN10139226_1.fastq /data/proctordm/cauris_genomes/fastq/SRA/SAMN10139226_2.fastq IL_B12938_Axilla_2019.sam
bash bwa.sh /data/proctordm/cauris_genomes/fastq/SRA/SAMN15689540_1.fastq /data/proctordm/cauris_genomes/fastq/SRA/SAMN15689540_2.fastq IL_CA01_Pleuralfluid_2021.sam
bash bwa.sh /data/proctordm/cauris_genomes/fastq/SRA/SAMN15689542_1.fastq /data/proctordm/cauris_genomes/fastq/SRA/SAMN15689542_2.fastq IL_CA03_urine_2021.sam
bash bwa.sh /data/proctordm/cauris_genomes/fastq/SRA/SAMN15689543_1.fastq /data/proctordm/cauris_genomes/fastq/SRA/SAMN15689543_2.fastq IL_CA04_respiratory_2021.sam
bash bwa.sh /data/proctordm/cauris_genomes/fastq/SRA/SAMN15689544_1.fastq /data/proctordm/cauris_genomes/fastq/SRA/SAMN15689544_2.fastq IL_CA05_blood_2021.sam
bash bwa.sh /data/proctordm/cauris_genomes/fastq/SRA/SAMN15689545_1.fastq /data/proctordm/cauris_genomes/fastq/SRA/SAMN15689545_2.fastq IL_CA06_blood_2021.sam
bash bwa.sh /data/proctordm/cauris_genomes/fastq/SRA/SAMN15689546_1.fastq /data/proctordm/cauris_genomes/fastq/SRA/SAMN15689546_2.fastq IL_CA07_wound_2021.sam
bash bwa.sh /data/proctordm/cauris_genomes/fastq/SRA/SAMN15689547_1.fastq /data/proctordm/cauris_genomes/fastq/SRA/SAMN15689547_2.fastq IL_CA08_Patientroom_2021.sam
bash bwa.sh /data/proctordm/cauris_genomes/fastq/SRA/SRR12363157_1.fastq /data/proctordm/cauris_genomes/fastq/SRA/SRR12363157_2.fastq IL_CA09_Patientroom_2021.sam
bash bwa.sh /data/proctordm/cauris_genomes/fastq/SRA/SAMN15689549_1.fastq /data/proctordm/cauris_genomes/fastq/SRA/SAMN15689549_2.fastq IL_CA10_Patientroom_2021.sam
bash bwa.sh /data/proctordm/cauris_genomes/fastq/SRA/SAMN15689550_1.fastq /data/proctordm/cauris_genomes/fastq/SRA/SAMN15689550_2.fastq IL_CA11_urine_2021.sam
bash bwa.sh /data/proctordm/cauris_genomes/fastq/SRA/SAMN15689552_1.fastq /data/proctordm/cauris_genomes/fastq/SRA/SAMN15689552_2.fastq IL_CA13_wound_2021.sam
bash bwa.sh /data/proctordm/cauris_genomes/fastq/SRA/SAMN15689553_1.fastq /data/proctordm/cauris_genomes/fastq/SRA/SAMN15689553_2.fastq IL_CA14_urine_2021.sam
bash bwa.sh /data/proctordm/cauris_genomes/fastq/SRA/SAMN15689554_1.fastq /data/proctordm/cauris_genomes/fastq/SRA/SAMN15689554_2.fastq IL_CA15_skin_2021.sam
bash bwa.sh /data/proctordm/cauris_genomes/fastq/SRA/SAMN15689555_1.fastq /data/proctordm/cauris_genomes/fastq/SRA/SAMN15689555_2.fastq IL_CA16_Nose_2021.sam
bash bwa.sh /data/proctordm/cauris_genomes/fastq/SRA/SAMN15689556_1.fastq /data/proctordm/cauris_genomes/fastq/SRA/SAMN15689556_2.fastq IL_CA17_Axilla_2021.sam
bash bwa.sh /data/proctordm/cauris_genomes/fastq/SRA/SAMN15689557_1.fastq /data/proctordm/cauris_genomes/fastq/SRA/SAMN15689557_2.fastq IL_CA18_Groin_2021.sam
bash bwa.sh /data/proctordm/cauris_genomes/fastq/SRA/SAMN15689558_1.fastq /data/proctordm/cauris_genomes/fastq/SRA/SAMN15689558_2.fastq IL_CA19_skin_2021.sam
bash bwa.sh /data/proctordm/cauris_genomes/fastq/SRA/SAMN15689559_1.fastq /data/proctordm/cauris_genomes/fastq/SRA/SAMN15689559_2.fastq IL_CA20_skin_2021.sam
bash bwa.sh /data/proctordm/cauris_genomes/fastq/SRA/SAMN15689560_1.fastq /data/proctordm/cauris_genomes/fastq/SRA/SAMN15689560_2.fastq IL_CA21_respiratory_2021.sam
bash bwa.sh /data/proctordm/cauris_genomes/fastq/SRA/SAMN15689561_1.fastq /data/proctordm/cauris_genomes/fastq/SRA/SAMN15689561_2.fastq IL_CA22_blood_2021.sam
bash bwa.sh /data/proctordm/cauris_genomes/fastq/SRA/SAMN15689562_1.fastq /data/proctordm/cauris_genomes/fastq/SRA/SAMN15689562_2.fastq IL_CA23_urine_2021.sam
bash bwa.sh /data/proctordm/cauris_genomes/fastq/SRA/SAMN15689563_1.fastq /data/proctordm/cauris_genomes/fastq/SRA/SAMN15689563_2.fastq IL_CA24_urine_2021.sam
bash bwa.sh /data/proctordm/cauris_genomes/fastq/SRA/SAMN15689564_1.fastq /data/proctordm/cauris_genomes/fastq/SRA/SAMN15689564_2.fastq IL_CA25_Synovialfluid_2021.sam
bash bwa.sh /data/Segrelab/data/all_genomes_reads/Caur_B15160_adqhu.BCC2418_1_20211028_182743_R1.fq.gz /data/Segrelab/data/all_genomes_reads/Caur_B15160_adqhu.BCC2418_1_20211028_182743_R2.fq.gz K_15160_Urine_2018.sam
bash bwa.sh /data/Segrelab/data/all_genomes_reads/Caur_B12189_adtwi.R1.fastq.gz /data/Segrelab/data/all_genomes_reads/Caur_B12189_adtwi.R2.fastq.gz K_B12189_Nares_2016.sam
bash bwa.sh /data/Segrelab/data/all_genomes_reads/Caur_B12510_adqhw.BCC2420_1_20211028_182743_R1.fq.gz /data/Segrelab/data/all_genomes_reads/Caur_B12510_adqhw.BCC2420_1_20211028_182743_R2.fq.gz K_B12510_AxGroin_2017.sam
bash bwa.sh /data/Segrelab/data/all_genomes_reads/Caur_B14707_adqht.BCC2417_1_20211028_182743_R1.fq.gz /data/Segrelab/data/all_genomes_reads/Caur_B14707_adqht.BCC2417_1_20211028_182743_R2.fq.gz K_B14707_AxGroin_2018.sam
bash bwa.sh /data/Segrelab/data/all_genomes_reads/Caur_B14714_adqhs.BCC2416_1_20211028_182743_R1.fq.gz /data/Segrelab/data/all_genomes_reads/Caur_B14714_adqhs.BCC2416_1_20211028_182743_R2.fq.gz K_B14714_AxGroin_2018.sam
bash bwa.sh /data/Segrelab/data/all_genomes_reads/Caur_B15161_adqhx.BCC2421_1_20211028_182743_R1.fq.gz /data/Segrelab/data/all_genomes_reads/Caur_B15161_adqhx.BCC2421_1_20211028_182743_R2.fq.gz K_B15161_Blood_2018.sam
bash bwa.sh /data/Segrelab/data/all_genomes_reads/Caur_B16276_adqhv.BCC2419_1_20211028_182743_R1.fq.gz /data/Segrelab/data/all_genomes_reads/Caur_B16276_adqhv.BCC2419_1_20211028_182743_R2.fq.gz K_B16276_Blood_2018.sam
bash bwa.sh /data/Segrelab/data/all_genomes_reads/Caur_B16278_adqhy.BCC2422_1_20211028_182743_R1.fq.gz /data/Segrelab/data/all_genomes_reads/Caur_B16278_adqhy.BCC2422_1_20211028_182743_R2.fq.gz K_B16278_Blood_2018.sam
bash bwa.sh /data/proctordm/cauris_genomes/fastq/SRA/SAMN10139374_1.fastq /data/proctordm/cauris_genomes/fastq/SRA/SAMN10139374_2.fastq MA_B12493_Bronchus_2019.sam
bash bwa.sh /data/Segrelab/data/all_genomes_reads/Caur.1.22.An_adtwq.R1.fastq.gz /data/Segrelab/data/all_genomes_reads/Caur.1.22.An_adtwq.R2.fastq.gz Sub14_An_1.sam
bash bwa.sh /data/Segrelab/data/all_genomes_reads/Caur.2.84.An_adtxf.R1.fastq.gz /data/Segrelab/data/all_genomes_reads/Caur.2.84.An_adtxf.R2.fastq.gz Sub14_An_2.sam
bash bwa.sh /data/Segrelab/data/all_genomes_reads/Caur.1.22.Fg_adtwr.R1.fastq.gz /data/Segrelab/data/all_genomes_reads/Caur.1.22.Fg_adtwr.R2.fastq.gz Sub14_Fg_1.sam
bash bwa.sh /data/Segrelab/data/all_genomes_reads/Caur.2.84.Fg_adtxg.R1.fastq.gz /data/Segrelab/data/all_genomes_reads/Caur.2.84.Fg_adtxg.R2.fastq.gz Sub14_Fg_2.sam
bash bwa.sh /data/Segrelab/data/all_genomes_reads/Caur.3.175.Ic_adtxt.R1.fastq.gz /data/Segrelab/data/all_genomes_reads/Caur.3.175.Ic_adtxt.R2.fastq.gz Sub14_Ic_3.sam
bash bwa.sh /data/Segrelab/data/all_genomes_reads/Caur.1.22.N_adtws.R1.fastq.gz /data/Segrelab/data/all_genomes_reads/Caur.1.22.N_adtws.R2.fastq.gz Sub14_N_1.sam
bash bwa.sh /data/Segrelab/data/all_genomes_reads/Caur.2.84.N_adtxh.R1.fastq.gz /data/Segrelab/data/all_genomes_reads/Caur.2.84.N_adtxh.R2.fastq.gz Sub14_N_2.sam
bash bwa.sh /data/Segrelab/data/all_genomes_reads/Caur.3.175.N_adtxu.R1.fastq.gz /data/Segrelab/data/all_genomes_reads/Caur.3.175.N_adtxu.R2.fastq.gz Sub14_N_3.sam
bash bwa.sh /data/Segrelab/data/all_genomes_reads/Caur.1.22.Tw_adtwt.R1.fastq.gz /data/Segrelab/data/all_genomes_reads/Caur.1.22.Tw_adtwt.R2.fastq.gz Sub14_Tw_1.sam
bash bwa.sh /data/Segrelab/data/all_genomes_reads/Caur.1.24.An_adtww.R1.fastq.gz /data/Segrelab/data/all_genomes_reads/Caur.1.24.An_adtww.R2.fastq.gz Sub15_An_1.sam
bash bwa.sh /data/Segrelab/data/all_genomes_reads/Caur.2.85.Fg_adtxi.R1.fastq.gz /data/Segrelab/data/all_genomes_reads/Caur.2.85.Fg_adtxi.R2.fastq.gz Sub15_Fg_2.sam
bash bwa.sh /data/Segrelab/data/all_genomes_reads/Caur.3.165.Fg_adtxr.R1.fastq.gz /data/Segrelab/data/all_genomes_reads/Caur.3.165.Fg_adtxr.R2.fastq.gz Sub15_Fg_3.sam
bash bwa.sh /data/Segrelab/data/all_genomes_reads/Caur.1.24.Ic_adtwx.R1.fastq.gz /data/Segrelab/data/all_genomes_reads/Caur.1.24.Ic_adtwx.R2.fastq.gz Sub15_Ic_1.sam
bash bwa.sh /data/Segrelab/data/all_genomes_reads/Caur.2.85.Tw_adtxj.R1.fastq.gz /data/Segrelab/data/all_genomes_reads/Caur.2.85.Tw_adtxj.R2.fastq.gz Sub15_Tw_2.sam
bash bwa.sh /data/Segrelab/data/all_genomes_reads/Caur.3.165.Tw_adtxs.R1.fastq.gz /data/Segrelab/data/all_genomes_reads/Caur.3.165.Tw_adtxs.R2.fastq.gz Sub15_Tw_3.sam
bash bwa.sh /data/Segrelab/data/all_genomes_reads/Caur.1.12.Ax_adtyf.R1.fastq.gz /data/Segrelab/data/all_genomes_reads/Caur.1.12.Ax_adtyf.R2.fastq.gz Sub2_Ax_1.sam
bash bwa.sh /data/Segrelab/data/all_genomes_reads/Caur.1.12.Fg_adtwj.R1.fastq.gz /data/Segrelab/data/all_genomes_reads/Caur.1.12.Fg_adtwj.R2.fastq.gz Sub2_Fg_1.sam
bash bwa.sh /data/Segrelab/data/all_genomes_reads/Caur.2.97.Fg_adtxk.R1.fastq.gz /data/Segrelab/data/all_genomes_reads/Caur.2.97.Fg_adtxk.R2.fastq.gz Sub2_Fg_2.sam
bash bwa.sh /data/Segrelab/data/all_genomes_reads/Caur.3.154.Fg_adtxp.R1.fastq.gz /data/Segrelab/data/all_genomes_reads/Caur.3.154.Fg_adtxp.R2.fastq.gz Sub2_Fg_3.sam
bash bwa.sh /data/Segrelab/data/all_genomes_reads/Caur.2.97.Ic_adtxl.R1.fastq.gz /data/Segrelab/data/all_genomes_reads/Caur.2.97.Ic_adtxl.R2.fastq.gz Sub2_Ic_2.sam
bash bwa.sh /data/Segrelab/data/all_genomes_reads/Caur.2.97.N_adtxm.R1.fastq.gz /data/Segrelab/data/all_genomes_reads/Caur.2.97.N_adtxm.R2.fastq.gz Sub2_N_2.sam
bash bwa.sh /data/Segrelab/data/all_genomes_reads/Caur.3.154.N_adtxq.R1.fastq.gz /data/Segrelab/data/all_genomes_reads/Caur.3.154.N_adtxq.R2.fastq.gz Sub2_N_3.sam
bash bwa.sh /data/Segrelab/data/all_genomes_reads/Caur.1.25.An_adtwy.R1.fastq.gz /data/Segrelab/data/all_genomes_reads/Caur.1.25.An_adtwy.R2.fastq.gz Sub23_An_1.sam
bash bwa.sh /data/Segrelab/data/all_genomes_reads/acrtx.AYD4730_1_20190730_152815_R1.fq.gz /data/Segrelab/data/all_genomes_reads/acrtx.AYD4730_1_20190730_152815_R2.fq.gz Sub23_Ea_2.sam
bash bwa.sh /data/Segrelab/data/all_genomes_reads/acuzs.AYD4744_1_20190730_152528_R1.fq.gz /data/Segrelab/data/all_genomes_reads/acuzs.AYD4744_1_20190730_152528_R2.fq.gz Sub23_Ea_3.sam
bash bwa.sh /data/Segrelab/data/all_genomes_reads/Caur.1.25.Fg_adtwz.R1.fastq.gz /data/Segrelab/data/all_genomes_reads/Caur.1.25.Fg_adtwz.R2.fastq.gz Sub23_Fg_1.sam
bash bwa.sh /data/Segrelab/data/all_genomes_reads/acuzt.AYD4745_1_20190730_153122_R1.fq.gz /data/Segrelab/data/all_genomes_reads/acuzt.AYD4745_1_20190730_153122_R2.fq.gz Sub23_Ic_3.sam
bash bwa.sh /data/Segrelab/data/all_genomes_reads/Caur.1.25.N_adtxa.R1.fastq.gz /data/Segrelab/data/all_genomes_reads/Caur.1.25.N_adtxa.R2.fastq.gz Sub23_N_1.sam
bash bwa.sh /data/Segrelab/data/all_genomes_reads/acuzg.AYD4732_1_20190730_155626_R1.fq.gz /data/Segrelab/data/all_genomes_reads/acuzg.AYD4732_1_20190730_155626_R2.fq.gz Sub23_N_2.sam
bash bwa.sh /data/Segrelab/data/all_genomes_reads/acuzu.AYD4746_1_20190730_153737_R1.fq.gz /data/Segrelab/data/all_genomes_reads/acuzu.AYD4746_1_20190730_153737_R2.fq.gz Sub23_N_3.sam
bash bwa.sh /data/Segrelab/data/all_genomes_reads/Caur.1.11.Fg_adtyd.R1.fastq.gz /data/Segrelab/data/all_genomes_reads/Caur.1.11.Fg_adtyd.R2.fastq.gz Sub28_Fg_1.sam
bash bwa.sh /data/Segrelab/data/all_genomes_reads/acuzk.AYD4736_1_20190730_165251_R1.fq.gz /data/Segrelab/data/all_genomes_reads/acuzk.AYD4736_1_20190730_165251_R2.fq.gz Sub28_Fg_3.sam
bash bwa.sh /data/Segrelab/data/all_genomes_reads/Caur.1.11.N_adtye.R1.fastq.gz /data/Segrelab/data/all_genomes_reads/Caur.1.11.N_adtye.R2.fastq.gz Sub28_N_1.sam
bash bwa.sh /data/Segrelab/data/all_genomes_reads/acrtr.AYD4724_1_20190731_042823_R1.fq.gz /data/Segrelab/data/all_genomes_reads/acrtr.AYD4724_1_20190731_042823_R2.fq.gz Sub28_N_2.sam
bash bwa.sh /data/Segrelab/data/all_genomes_reads/acuzl.AYD4737_1_20190730_165948_R1.fq.gz /data/Segrelab/data/all_genomes_reads/acuzl.AYD4737_1_20190730_165948_R2.fq.gz Sub28_N_3.sam
bash bwa.sh /data/Segrelab/data/all_genomes_reads/acuzv.AYD4747_1_20190730_154305_R1.fq.gz /data/Segrelab/data/all_genomes_reads/acuzv.AYD4747_1_20190730_154305_R2.fq.gz Sub35_An_3.sam
bash bwa.sh /data/Segrelab/data/all_genomes_reads/acuzi.AYD4734_1_20190730_161754_R1.fq.gz /data/Segrelab/data/all_genomes_reads/acuzi.AYD4734_1_20190730_161754_R2.fq.gz Sub35_Ic_2.sam
bash bwa.sh /data/Segrelab/data/all_genomes_reads/acuzw.AYD4748_1_20190730_154847_R1.fq.gz /data/Segrelab/data/all_genomes_reads/acuzw.AYD4748_1_20190730_154847_R2.fq.gz Sub35_Ic_3.sam
bash bwa.sh /data/Segrelab/data/all_genomes_reads/acuzj.AYD4735_1_20190730_164546_R1.fq.gz /data/Segrelab/data/all_genomes_reads/acuzj.AYD4735_1_20190730_164546_R2.fq.gz Sub35_N_2.sam
bash bwa.sh /data/Segrelab/data/all_genomes_reads/acuzx.AYD4749_1_20190730_155601_R1.fq.gz /data/Segrelab/data/all_genomes_reads/acuzx.AYD4749_1_20190730_155601_R2.fq.gz Sub35_N_3.sam
bash bwa.sh /data/Segrelab/data/all_genomes_reads/Caur.1.13.Ic_adtwk.R1.fastq.gz /data/Segrelab/data/all_genomes_reads/Caur.1.13.Ic_adtwk.R2.fastq.gz Sub4_Ic_1.sam
bash bwa.sh /data/Segrelab/data/all_genomes_reads/Caur.3.185.Ic_adtxw.R1.fastq.gz /data/Segrelab/data/all_genomes_reads/Caur.3.185.Ic_adtxw.R2.fastq.gz Sub4_Ic_3.sam
bash bwa.sh /data/Segrelab/data/all_genomes_reads/Caur.1.13.N_adtyg.R1.fastq.gz /data/Segrelab/data/all_genomes_reads/Caur.1.13.N_adtyg.R2.fastq.gz Sub4_N_1.sam
bash bwa.sh /data/Segrelab/data/all_genomes_reads/Caur.2.98.N_adtxn.R1.fastq.gz /data/Segrelab/data/all_genomes_reads/Caur.2.98.N_adtxn.R2.fastq.gz Sub4_N_2.sam
bash bwa.sh /data/Segrelab/data/all_genomes_reads/Caur.3.185.N_adtxx.R1.fastq.gz /data/Segrelab/data/all_genomes_reads/Caur.3.185.N_adtxx.R2.fastq.gz Sub4_N_3.sam
bash bwa.sh /data/Segrelab/data/all_genomes_reads/Caur.2.98.Tw_adtxo.R1.fastq.gz /data/Segrelab/data/all_genomes_reads/Caur.2.98.Tw_adtxo.R2.fastq.gz Sub4_Tw_2.sam
bash bwa.sh /data/Segrelab/data/all_genomes_reads/Caur.1.10.An_adtxb.R1.fastq.gz /data/Segrelab/data/all_genomes_reads/Caur.1.10.An_adtxb.R2.fastq.gz Sub46_An_1.sam
bash bwa.sh /data/Segrelab/data/all_genomes_reads/Caur.1.10.Fg_adtyb.R1.fastq.gz /data/Segrelab/data/all_genomes_reads/Caur.1.10.Fg_adtyb.R2.fastq.gz Sub46_Fg_1.sam
bash bwa.sh /data/Segrelab/data/all_genomes_reads/acrts.AYD4725_1_20190731_043616_R1.fq.gz /data/Segrelab/data/all_genomes_reads/acrts.AYD4725_1_20190731_043616_R2.fq.gz Sub46_Fg_2.sam
bash bwa.sh /data/Segrelab/data/all_genomes_reads/acuzm.AYD4738_1_20190730_170811_R1.fq.gz /data/Segrelab/data/all_genomes_reads/acuzm.AYD4738_1_20190730_170811_R2.fq.gz Sub46_Fg_3.sam
bash bwa.sh /data/Segrelab/data/all_genomes_reads/Caur.1.10.N_adtyc.R1.fastq.gz /data/Segrelab/data/all_genomes_reads/Caur.1.10.N_adtyc.R2.fastq.gz Sub46_N_1.sam
bash bwa.sh /data/Segrelab/data/all_genomes_reads/acrtt.AYD4726_1_20190731_044304_R1.fq.gz /data/Segrelab/data/all_genomes_reads/acrtt.AYD4726_1_20190731_044304_R2.fq.gz Sub46_N_2.sam
bash bwa.sh /data/Segrelab/data/all_genomes_reads/acuzn.AYD4739_1_20190730_145801_R1.fq.gz /data/Segrelab/data/all_genomes_reads/acuzn.AYD4739_1_20190730_145801_R2.fq.gz Sub46_N_3.sam
bash bwa.sh /data/Segrelab/data/all_genomes_reads/acrtu.AYD4727_1_20190731_045125_R1.fq.gz /data/Segrelab/data/all_genomes_reads/acrtu.AYD4727_1_20190731_045125_R2.fq.gz Sub46_Ne_2.sam
bash bwa.sh /data/Segrelab/data/all_genomes_reads/acuzo.AYD4740_1_20190730_150404_R1.fq.gz /data/Segrelab/data/all_genomes_reads/acuzo.AYD4740_1_20190730_150404_R2.fq.gz Sub46_Ne_3.sam
bash bwa.sh /data/Segrelab/data/all_genomes_reads/Caur.1.14.An_adtwl.R1.fastq.gz /data/Segrelab/data/all_genomes_reads/Caur.1.14.An_adtwl.R2.fastq.gz Sub48_An_1.sam
bash bwa.sh /data/Segrelab/data/all_genomes_reads/Caur.3.186.An_adtxy.R1.fastq.gz /data/Segrelab/data/all_genomes_reads/Caur.3.186.An_adtxy.R2.fastq.gz Sub48_An_3.sam
bash bwa.sh /data/Segrelab/data/all_genomes_reads/Caur.1.14.Fg_adtwm.R1.fastq.gz /data/Segrelab/data/all_genomes_reads/Caur.1.14.Fg_adtwm.R2.fastq.gz Sub48_Fg_1.sam
bash bwa.sh /data/Segrelab/data/all_genomes_reads/Caur.3.186.Fg_adtxz.R1.fastq.gz /data/Segrelab/data/all_genomes_reads/Caur.3.186.Fg_adtxz.R2.fastq.gz Sub48_Fg_3.sam
bash bwa.sh /data/Segrelab/data/all_genomes_reads/Caur.1.14.Ic_adtwn.R1.fastq.gz /data/Segrelab/data/all_genomes_reads/Caur.1.14.Ic_adtwn.R2.fastq.gz Sub48_Ic_1.sam
bash bwa.sh /data/Segrelab/data/all_genomes_reads/Caur.1.14.N_adtwo.R1.fastq.gz /data/Segrelab/data/all_genomes_reads/Caur.1.14.N_adtwo.R2.fastq.gz Sub48_N_1.sam
bash bwa.sh /data/Segrelab/data/all_genomes_reads/Caur.3.186.N_adtya.R1.fastq.gz /data/Segrelab/data/all_genomes_reads/Caur.3.186.N_adtya.R2.fastq.gz Sub48_N_3.sam
bash bwa.sh /data/Segrelab/data/all_genomes_reads/Caur.1.14.Tw_adtwp.R1.fastq.gz /data/Segrelab/data/all_genomes_reads/Caur.1.14.Tw_adtwp.R2.fastq.gz Sub48_Tw_1.sam
bash bwa.sh /data/Segrelab/data/all_genomes_reads/Caur.1.23.An_adtwu.R1.fastq.gz /data/Segrelab/data/all_genomes_reads/Caur.1.23.An_adtwu.R2.fastq.gz Sub5_An_1.sam
bash bwa.sh /data/Segrelab/data/all_genomes_reads/Caur.2.82.Ea_adtxe.R1.fastq.gz /data/Segrelab/data/all_genomes_reads/Caur.2.82.Ea_adtxe.R2.fastq.gz Sub5_Ea_2.sam
bash bwa.sh /data/Segrelab/data/all_genomes_reads/Caur.1.23.N_adtwv.R1.fastq.gz /data/Segrelab/data/all_genomes_reads/Caur.1.23.N_adtwv.R2.fastq.gz Sub5_N_1.sam
bash bwa.sh /data/Segrelab/data/all_genomes_reads/Caur.3.176.N_adtxv.R1.fastq.gz /data/Segrelab/data/all_genomes_reads/Caur.3.176.N_adtxv.R2.fastq.gz Sub5_N_3.sam
bash bwa.sh /data/Segrelab/data/all_genomes_reads/Caur.1.53.Fg_adtxc.R1.fastq.gz /data/Segrelab/data/all_genomes_reads/Caur.1.53.Fg_adtxc.R2.fastq.gz Sub53_Fg_1.sam
bash bwa.sh /data/Segrelab/data/all_genomes_reads/acuzp.AYD4741_1_20190730_150938_R1.fq.gz /data/Segrelab/data/all_genomes_reads/acuzp.AYD4741_1_20190730_150938_R2.fq.gz Sub53_Fg_3.sam
bash bwa.sh /data/Segrelab/data/all_genomes_reads/Caur.1.53.N_adtxd.R1.fastq.gz /data/Segrelab/data/all_genomes_reads/Caur.1.53.N_adtxd.R2.fastq.gz Sub53_N_1.sam
bash bwa.sh /data/Segrelab/data/all_genomes_reads/acrtv.AYD4728_1_20190730_145728_R1.fq.gz /data/Segrelab/data/all_genomes_reads/acrtv.AYD4728_1_20190730_145728_R2.fq.gz Sub53_N_2.sam
bash bwa.sh /data/Segrelab/data/all_genomes_reads/acuzq.AYD4742_1_20190730_151511_R1.fq.gz /data/Segrelab/data/all_genomes_reads/acuzq.AYD4742_1_20190730_151511_R2.fq.gz Sub53_N_3.sam
bash bwa.sh /data/Segrelab/data/all_genomes_reads/acrtw.AYD4729_1_20190730_151525_R1.fq.gz /data/Segrelab/data/all_genomes_reads/acrtw.AYD4729_1_20190730_151525_R2.fq.gz Sub53_Tw_2.sam
bash bwa.sh /data/Segrelab/data/all_genomes_reads/acuzr.AYD4743_1_20190730_152013_R1.fq.gz /data/Segrelab/data/all_genomes_reads/acuzr.AYD4743_1_20190730_152013_R2.fq.gz Sub53_Tw_3.sam
bash bwa.sh /data/Segrelab/data/all_genomes_reads/Caur_B11889_adqhz.BCC2423_1_20211028_182743_R1.fq.gz /data/Segrelab/data/all_genomes_reads/Caur_B11889_adqhz.BCC2423_1_20211028_182743_R2.fq.gz UC_B11889_Env_2016.sam
bash bwa.sh /data/Segrelab/data/all_genomes_reads/Caur_B11891_adqia.BCC2424_1_20211028_182743_R1.fq.gz /data/Segrelab/data/all_genomes_reads/Caur_B11891_adqia.BCC2424_1_20211028_182743_R2.fq.gz UC_B11891_Env_2016.sam
swarm -f bwa.swarm --job-name bwa -t 20 -g 128 --module bwa --time 4:00:00 --logdir bwa_swarm_out 
mv bwa_swarm_out ../scripts

cat <<EOF >gettotals.sh
#!/usr/bin/bash
bam_files=\$1 
flagstat=\$2
module load samtools
samtools flagstat \${bam_files} > \${flagstat}
EOF

###generate the swarm file
bam_files=($(ls  *.sam ))
flagstat=($(echo "${bam_files[*]}" | sed 's/.sam/_flagstat.txt/g'))

total=${#bam_files[*]}
for ((i=0; i<=$(($total -1)); i++))
do
	echo "bash gettotals.sh ${bam_files[i]} ${flagstat[i]}"
done > totals.swarm
swarm -f totals.swarm --job-name totals -t 1 -g 128 --time 2:00:00 --logdir /data/proctordm/cauris_genomes/scripts/mapping_totals



######################### Run Picard to clean up bam filess
#split thee filees into subdirectories of size 5
n=0; for f in *; do d="subdir$((n++ / 2))"; mkdir -p "$d"; mv -- "$f" "$d/$f"; done

for d in */; do cp picard.sh "$d"; done


###### Execute your script in on all  sudirectories   
mydir=$PWD
for d in */; 
    do
        cd $d
        sbatch --mem=98g --cpus-per-task=16 --time 5:00:00 --gres=lscratch:400 picard.sh 
        cd $mydir;
    done



####run picardv5.0.sh
#!/usr/bin/bash
#read this tutorial on population genomics; https://informatics.fas.harvard.edu/whole-genome-resquencing-for-population-genomics-fastq-to-vcf.html
#submit with: sbatch --cpus-per-task=16 --mem=128g --gres=lscratch:200 --time 24:00:00 picardv5.0.sh
#relative to v3.0 fixed the name for 'derep_bam_names_sorted' which is causing the output to fail, I think
#relative to v5.0 I delete files here and don't validate sam since i know we're missing read groups
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
        java -Xmx96g -jar $PICARDJARPATH/picard.jar CollectAlignmentSummaryMetrics I=${sorted_bam_names[i]} R=//data/proctordm/cauris_genomes/reference/Caur.B11891_adqia.Nano.pilon.fasta METRIC_ACCUMULATION_LEVEL=SAMPLE METRIC_ACCUMULATION_LEVEL=READ_GROUP O=${picard_alignment_metrics[i]}
        java -Xmx96g -jar $PICARDJARPATH/picard.jar MarkDuplicates I=${sorted_bam_names[i]} O=${derep_bam_names[i]} METRICS_FILE=${picard_derep_text[i]} VALIDATION_STRINGENCY=SILENT 
        java -Xmx96g -jar $PICARDJARPATH/picard.jar SortSam I=${derep_bam_names[i]} O=${derep_bam_names_sorted[i]} SORT_ORDER=coordinate VALIDATION_STRINGENCY=SILENT CREATE_INDEX=false  
        rm ${sorted_bam_names[i]}
        rm ${derep_bam_names[i]}
done

######copy the picard output to the main directory
mv subdir*/* .



####run picard to add read groups
cat << EOF > addRG.sh
#!/usr/bin/bash
module load picard
module load samtools
BAMIN=\$1 
BAMOUT=\$2
RG=\$3

java -Xmx96g -jar \$PICARDJARPATH/picard.jar AddOrReplaceReadGroups \
       I=\${BAMIN}\
       O=\${BAMOUT}\
       RGID=\${RG}\
       RGLB="none"\
       RGPL="none"\
       RGPU="none"\
       RGSM=\${RG}
EOF


######################### Run Picard to clean up bam filess
nano addrg.swarm
bash addRG.sh FL_B12847_urine_2019_derep_sorted.bam rg_out/FL_B12847_urine_2019.bam FL_B12847_urine_2019
bash addRG.sh FL_B12847_urine_2021_derep_sorted.bam rg_out/FL_B12847_urine_2021.bam FL_B12847_urine_2021
bash addRG.sh FL_B16329_blood_2020_derep_sorted.bam rg_out/FL_B16329_blood_2020.bam FL_B16329_blood_2020
bash addRG.sh IL_B11842_urine_2019_derep_sorted.bam rg_out/IL_B11842_urine_2019.bam IL_B11842_urine_2019
bash addRG.sh IL_B11843_blood_2019_derep_sorted.bam rg_out/IL_B11843_blood_2019.bam IL_B11843_blood_2019
bash addRG.sh IL_B11880_nares_2019_derep_sorted.bam rg_out/IL_B11880_nares_2019.bam IL_B11880_nares_2019
bash addRG.sh IL_B11882_nares_2019_derep_sorted.bam rg_out/IL_B11882_nares_2019.bam IL_B11882_nares_2019
bash addRG.sh IL_B11883_Rectal_2019_derep_sorted.bam rg_out/IL_B11883_Rectal_2019.bam IL_B11883_Rectal_2019
bash addRG.sh IL_B11884_Vaginal_2019_derep_sorted.bam rg_out/IL_B11884_Vaginal_2019.bam IL_B11884_Vaginal_2019
bash addRG.sh IL_B11885_Axilla_2019_derep_sorted.bam rg_out/IL_B11885_Axilla_2019.bam IL_B11885_Axilla_2019
bash addRG.sh IL_B11886_Groin_2019_derep_sorted.bam rg_out/IL_B11886_Groin_2019.bam IL_B11886_Groin_2019
bash addRG.sh IL_B12028_blood_2019_derep_sorted.bam rg_out/IL_B12028_blood_2019.bam IL_B12028_blood_2019
bash addRG.sh IL_B12029_blood_2019_derep_sorted.bam rg_out/IL_B12029_blood_2019.bam IL_B12029_blood_2019
bash addRG.sh IL_B12030_blood_2019_derep_sorted.bam rg_out/IL_B12030_blood_2019.bam IL_B12030_blood_2019
bash addRG.sh IL_B12031_blood_2019_derep_sorted.bam rg_out/IL_B12031_blood_2019.bam IL_B12031_blood_2019
bash addRG.sh IL_B12032_blood_2019_derep_sorted.bam rg_out/IL_B12032_blood_2019.bam IL_B12032_blood_2019
bash addRG.sh IL_B12033_blood_2019_derep_sorted.bam rg_out/IL_B12033_blood_2019.bam IL_B12033_blood_2019
bash addRG.sh IL_B12034_blood_2019_derep_sorted.bam rg_out/IL_B12034_blood_2019.bam IL_B12034_blood_2019
bash addRG.sh IL_B12035_blood_2019_derep_sorted.bam rg_out/IL_B12035_blood_2019.bam IL_B12035_blood_2019
bash addRG.sh IL_B12036_blood_2019_derep_sorted.bam rg_out/IL_B12036_blood_2019.bam IL_B12036_blood_2019
bash addRG.sh IL_B12046_AxillaGroin_2019_derep_sorted.bam rg_out/IL_B12046_AxillaGroin_2019.bam IL_B12046_AxillaGroin_2019
bash addRG.sh IL_B12077_feces_2019_derep_sorted.bam rg_out/IL_B12077_feces_2019.bam IL_B12077_feces_2019
bash addRG.sh IL_B12189_nares_2019_derep_sorted.bam rg_out/IL_B12189_nares_2019.bam IL_B12189_nares_2019
bash addRG.sh IL_B12190_Axilla_2019_derep_sorted.bam rg_out/IL_B12190_Axilla_2019.bam IL_B12190_Axilla_2019
bash addRG.sh IL_B12229_kidney_2019_derep_sorted.bam rg_out/IL_B12229_kidney_2019.bam IL_B12229_kidney_2019
bash addRG.sh IL_B12378_AxillaGroin_2019_derep_sorted.bam rg_out/IL_B12378_AxillaGroin_2019.bam IL_B12378_AxillaGroin_2019
bash addRG.sh IL_B12380_AxillaGroin_2019_derep_sorted.bam rg_out/IL_B12380_AxillaGroin_2019.bam IL_B12380_AxillaGroin_2019
bash addRG.sh IL_B12388_AxillaGroin_2019_derep_sorted.bam rg_out/IL_B12388_AxillaGroin_2019.bam IL_B12388_AxillaGroin_2019
bash addRG.sh IL_B12406_urine_2019_derep_sorted.bam rg_out/IL_B12406_urine_2019.bam IL_B12406_urine_2019
bash addRG.sh IL_B12510_AxillaGroin_2019_derep_sorted.bam rg_out/IL_B12510_AxillaGroin_2019.bam IL_B12510_AxillaGroin_2019
bash addRG.sh IL_B12579_Axilla_2019_derep_sorted.bam rg_out/IL_B12579_Axilla_2019.bam IL_B12579_Axilla_2019
bash addRG.sh IL_B12937_Axilla_2019_derep_sorted.bam rg_out/IL_B12937_Axilla_2019.bam IL_B12937_Axilla_2019
bash addRG.sh IL_B12938_Axilla_2019_derep_sorted.bam rg_out/IL_B12938_Axilla_2019.bam IL_B12938_Axilla_2019
bash addRG.sh IL_CA01_Pleuralfluid_2021_derep_sorted.bam rg_out/IL_CA01_Pleuralfluid_2021.bam IL_CA01_Pleuralfluid_2021
bash addRG.sh IL_CA03_urine_2021_derep_sorted.bam rg_out/IL_CA03_urine_2021.bam IL_CA03_urine_2021
bash addRG.sh IL_CA04_respiratory_2021_derep_sorted.bam rg_out/IL_CA04_respiratory_2021.bam IL_CA04_respiratory_2021
bash addRG.sh IL_CA05_blood_2021_derep_sorted.bam rg_out/IL_CA05_blood_2021.bam IL_CA05_blood_2021
bash addRG.sh IL_CA06_blood_2021_derep_sorted.bam rg_out/IL_CA06_blood_2021.bam IL_CA06_blood_2021
bash addRG.sh IL_CA07_wound_2021_derep_sorted.bam rg_out/IL_CA07_wound_2021.bam IL_CA07_wound_2021
bash addRG.sh IL_CA08_Patientroom_2021_derep_sorted.bam rg_out/IL_CA08_Patientroom_2021.bam IL_CA08_Patientroom_2021
bash addRG.sh IL_CA09_Patientroom_2021_derep_sorted.bam rg_out/IL_CA09_Patientroom_2021.bam IL_CA09_Patientroom_2021
bash addRG.sh IL_CA10_Patientroom_2021_derep_sorted.bam rg_out/IL_CA10_Patientroom_2021.bam IL_CA10_Patientroom_2021
bash addRG.sh IL_CA11_urine_2021_derep_sorted.bam rg_out/IL_CA11_urine_2021.bam IL_CA11_urine_2021
bash addRG.sh IL_CA13_wound_2021_derep_sorted.bam rg_out/IL_CA13_wound_2021.bam IL_CA13_wound_2021
bash addRG.sh IL_CA14_urine_2021_derep_sorted.bam rg_out/IL_CA14_urine_2021.bam IL_CA14_urine_2021
bash addRG.sh IL_CA15_skin_2021_derep_sorted.bam rg_out/IL_CA15_skin_2021.bam IL_CA15_skin_2021
bash addRG.sh IL_CA16_Nose_2021_derep_sorted.bam rg_out/IL_CA16_Nose_2021.bam IL_CA16_Nose_2021
bash addRG.sh IL_CA17_Axilla_2021_derep_sorted.bam rg_out/IL_CA17_Axilla_2021.bam IL_CA17_Axilla_2021
bash addRG.sh IL_CA18_Groin_2021_derep_sorted.bam rg_out/IL_CA18_Groin_2021.bam IL_CA18_Groin_2021
bash addRG.sh IL_CA19_skin_2021_derep_sorted.bam rg_out/IL_CA19_skin_2021.bam IL_CA19_skin_2021
bash addRG.sh IL_CA20_skin_2021_derep_sorted.bam rg_out/IL_CA20_skin_2021.bam IL_CA20_skin_2021
bash addRG.sh IL_CA21_respiratory_2021_derep_sorted.bam rg_out/IL_CA21_respiratory_2021.bam IL_CA21_respiratory_2021
bash addRG.sh IL_CA22_blood_2021_derep_sorted.bam rg_out/IL_CA22_blood_2021.bam IL_CA22_blood_2021
bash addRG.sh IL_CA23_urine_2021_derep_sorted.bam rg_out/IL_CA23_urine_2021.bam IL_CA23_urine_2021
bash addRG.sh IL_CA24_urine_2021_derep_sorted.bam rg_out/IL_CA24_urine_2021.bam IL_CA24_urine_2021
bash addRG.sh IL_CA25_Synovialfluid_2021_derep_sorted.bam rg_out/IL_CA25_Synovialfluid_2021.bam IL_CA25_Synovialfluid_2021
bash addRG.sh K_15160_Urine_2018_derep_sorted.bam rg_out/K_15160_Urine_2018.bam K_15160_Urine_2018
bash addRG.sh K_B12189_Nares_2016_derep_sorted.bam rg_out/K_B12189_Nares_2016.bam K_B12189_Nares_2016
bash addRG.sh K_B12510_AxGroin_2017_derep_sorted.bam rg_out/K_B12510_AxGroin_2017.bam K_B12510_AxGroin_2017
bash addRG.sh K_B14707_AxGroin_2018_derep_sorted.bam rg_out/K_B14707_AxGroin_2018.bam K_B14707_AxGroin_2018
bash addRG.sh K_B14714_AxGroin_2018_derep_sorted.bam rg_out/K_B14714_AxGroin_2018.bam K_B14714_AxGroin_2018
bash addRG.sh K_B15161_Blood_2018_derep_sorted.bam rg_out/K_B15161_Blood_2018.bam K_B15161_Blood_2018
bash addRG.sh K_B16276_Blood_2018_derep_sorted.bam rg_out/K_B16276_Blood_2018.bam K_B16276_Blood_2018
bash addRG.sh K_B16278_Blood_2018_derep_sorted.bam rg_out/K_B16278_Blood_2018.bam K_B16278_Blood_2018
bash addRG.sh MA_B12493_Bronchus_2019_derep_sorted.bam rg_out/MA_B12493_Bronchus_2019.bam MA_B12493_Bronchus_2019
bash addRG.sh Sub14_An_1_derep_sorted.bam rg_out/Sub14_An_1.bam Sub14_An_1
bash addRG.sh Sub14_An_2_derep_sorted.bam rg_out/Sub14_An_2.bam Sub14_An_2
bash addRG.sh Sub14_Fg_1_derep_sorted.bam rg_out/Sub14_Fg_1.bam Sub14_Fg_1
bash addRG.sh Sub14_Fg_2_derep_sorted.bam rg_out/Sub14_Fg_2.bam Sub14_Fg_2
bash addRG.sh Sub14_Ic_3_derep_sorted.bam rg_out/Sub14_Ic_3.bam Sub14_Ic_3
bash addRG.sh Sub14_N_1_derep_sorted.bam rg_out/Sub14_N_1.bam Sub14_N_1
bash addRG.sh Sub14_N_2_derep_sorted.bam rg_out/Sub14_N_2.bam Sub14_N_2
bash addRG.sh Sub14_N_3_derep_sorted.bam rg_out/Sub14_N_3.bam Sub14_N_3
bash addRG.sh Sub14_Tw_1_derep_sorted.bam rg_out/Sub14_Tw_1.bam Sub14_Tw_1
bash addRG.sh Sub15_An_1_derep_sorted.bam rg_out/Sub15_An_1.bam Sub15_An_1
bash addRG.sh Sub15_Fg_2_derep_sorted.bam rg_out/Sub15_Fg_2.bam Sub15_Fg_2
bash addRG.sh Sub15_Fg_3_derep_sorted.bam rg_out/Sub15_Fg_3.bam Sub15_Fg_3
bash addRG.sh Sub15_Ic_1_derep_sorted.bam rg_out/Sub15_Ic_1.bam Sub15_Ic_1
bash addRG.sh Sub15_Tw_2_derep_sorted.bam rg_out/Sub15_Tw_2.bam Sub15_Tw_2
bash addRG.sh Sub15_Tw_3_derep_sorted.bam rg_out/Sub15_Tw_3.bam Sub15_Tw_3
bash addRG.sh Sub2_Ax_1_derep_sorted.bam rg_out/Sub2_Ax_1.bam Sub2_Ax_1
bash addRG.sh Sub2_Fg_1_derep_sorted.bam rg_out/Sub2_Fg_1.bam Sub2_Fg_1
bash addRG.sh Sub2_Fg_2_derep_sorted.bam rg_out/Sub2_Fg_2.bam Sub2_Fg_2
bash addRG.sh Sub2_Fg_3_derep_sorted.bam rg_out/Sub2_Fg_3.bam Sub2_Fg_3
bash addRG.sh Sub2_Ic_2_derep_sorted.bam rg_out/Sub2_Ic_2.bam Sub2_Ic_2
bash addRG.sh Sub2_N_2_derep_sorted.bam rg_out/Sub2_N_2.bam Sub2_N_2
bash addRG.sh Sub2_N_3_derep_sorted.bam rg_out/Sub2_N_3.bam Sub2_N_3
bash addRG.sh Sub23_An_1_derep_sorted.bam rg_out/Sub23_An_1.bam Sub23_An_1
bash addRG.sh Sub23_Ea_2_derep_sorted.bam rg_out/Sub23_Ea_2.bam Sub23_Ea_2
bash addRG.sh Sub23_Ea_3_derep_sorted.bam rg_out/Sub23_Ea_3.bam Sub23_Ea_3
bash addRG.sh Sub23_Fg_1_derep_sorted.bam rg_out/Sub23_Fg_1.bam Sub23_Fg_1
bash addRG.sh Sub23_Ic_3_derep_sorted.bam rg_out/Sub23_Ic_3.bam Sub23_Ic_3
bash addRG.sh Sub23_N_1_derep_sorted.bam rg_out/Sub23_N_1.bam Sub23_N_1
bash addRG.sh Sub23_N_2_derep_sorted.bam rg_out/Sub23_N_2.bam Sub23_N_2
bash addRG.sh Sub23_N_3_derep_sorted.bam rg_out/Sub23_N_3.bam Sub23_N_3
bash addRG.sh Sub28_Fg_1_derep_sorted.bam rg_out/Sub28_Fg_1.bam Sub28_Fg_1
bash addRG.sh Sub28_Fg_3_derep_sorted.bam rg_out/Sub28_Fg_3.bam Sub28_Fg_3
bash addRG.sh Sub28_N_1_derep_sorted.bam rg_out/Sub28_N_1.bam Sub28_N_1
bash addRG.sh Sub28_N_2_derep_sorted.bam rg_out/Sub28_N_2.bam Sub28_N_2
bash addRG.sh Sub28_N_3_derep_sorted.bam rg_out/Sub28_N_3.bam Sub28_N_3
bash addRG.sh Sub35_An_3_derep_sorted.bam rg_out/Sub35_An_3.bam Sub35_An_3
bash addRG.sh Sub35_Ic_2_derep_sorted.bam rg_out/Sub35_Ic_2.bam Sub35_Ic_2
bash addRG.sh Sub35_Ic_3_derep_sorted.bam rg_out/Sub35_Ic_3.bam Sub35_Ic_3
bash addRG.sh Sub35_N_2_derep_sorted.bam rg_out/Sub35_N_2.bam Sub35_N_2
bash addRG.sh Sub35_N_3_derep_sorted.bam rg_out/Sub35_N_3.bam Sub35_N_3
bash addRG.sh Sub4_Ic_1_derep_sorted.bam rg_out/Sub4_Ic_1.bam Sub4_Ic_1
bash addRG.sh Sub4_Ic_3_derep_sorted.bam rg_out/Sub4_Ic_3.bam Sub4_Ic_3
bash addRG.sh Sub4_N_1_derep_sorted.bam rg_out/Sub4_N_1.bam Sub4_N_1
bash addRG.sh Sub4_N_2_derep_sorted.bam rg_out/Sub4_N_2.bam Sub4_N_2
bash addRG.sh Sub4_N_3_derep_sorted.bam rg_out/Sub4_N_3.bam Sub4_N_3
bash addRG.sh Sub4_Tw_2_derep_sorted.bam rg_out/Sub4_Tw_2.bam Sub4_Tw_2
bash addRG.sh Sub46_An_1_derep_sorted.bam rg_out/Sub46_An_1.bam Sub46_An_1
bash addRG.sh Sub46_Fg_1_derep_sorted.bam rg_out/Sub46_Fg_1.bam Sub46_Fg_1
bash addRG.sh Sub46_Fg_2_derep_sorted.bam rg_out/Sub46_Fg_2.bam Sub46_Fg_2
bash addRG.sh Sub46_Fg_3_derep_sorted.bam rg_out/Sub46_Fg_3.bam Sub46_Fg_3
bash addRG.sh Sub46_N_1_derep_sorted.bam rg_out/Sub46_N_1.bam Sub46_N_1
bash addRG.sh Sub46_N_2_derep_sorted.bam rg_out/Sub46_N_2.bam Sub46_N_2
bash addRG.sh Sub46_N_3_derep_sorted.bam rg_out/Sub46_N_3.bam Sub46_N_3
bash addRG.sh Sub46_Ne_2_derep_sorted.bam rg_out/Sub46_Ne_2.bam Sub46_Ne_2
bash addRG.sh Sub46_Ne_3_derep_sorted.bam rg_out/Sub46_Ne_3.bam Sub46_Ne_3
bash addRG.sh Sub48_An_1_derep_sorted.bam rg_out/Sub48_An_1.bam Sub48_An_1
bash addRG.sh Sub48_An_3_derep_sorted.bam rg_out/Sub48_An_3.bam Sub48_An_3
bash addRG.sh Sub48_Fg_1_derep_sorted.bam rg_out/Sub48_Fg_1.bam Sub48_Fg_1
bash addRG.sh Sub48_Fg_3_derep_sorted.bam rg_out/Sub48_Fg_3.bam Sub48_Fg_3
bash addRG.sh Sub48_Ic_1_derep_sorted.bam rg_out/Sub48_Ic_1.bam Sub48_Ic_1
bash addRG.sh Sub48_N_1_derep_sorted.bam rg_out/Sub48_N_1.bam Sub48_N_1
bash addRG.sh Sub48_N_3_derep_sorted.bam rg_out/Sub48_N_3.bam Sub48_N_3
bash addRG.sh Sub48_Tw_1_derep_sorted.bam rg_out/Sub48_Tw_1.bam Sub48_Tw_1
bash addRG.sh Sub5_An_1_derep_sorted.bam rg_out/Sub5_An_1.bam Sub5_An_1
bash addRG.sh Sub5_Ea_2_derep_sorted.bam rg_out/Sub5_Ea_2.bam Sub5_Ea_2
bash addRG.sh Sub5_N_1_derep_sorted.bam rg_out/Sub5_N_1.bam Sub5_N_1
bash addRG.sh Sub5_N_3_derep_sorted.bam rg_out/Sub5_N_3.bam Sub5_N_3
bash addRG.sh Sub53_Fg_1_derep_sorted.bam rg_out/Sub53_Fg_1.bam Sub53_Fg_1
bash addRG.sh Sub53_Fg_3_derep_sorted.bam rg_out/Sub53_Fg_3.bam Sub53_Fg_3
bash addRG.sh Sub53_N_1_derep_sorted.bam rg_out/Sub53_N_1.bam Sub53_N_1
bash addRG.sh Sub53_N_2_derep_sorted.bam rg_out/Sub53_N_2.bam Sub53_N_2
bash addRG.sh Sub53_N_3_derep_sorted.bam rg_out/Sub53_N_3.bam Sub53_N_3
bash addRG.sh Sub53_Tw_2_derep_sorted.bam rg_out/Sub53_Tw_2.bam Sub53_Tw_2
bash addRG.sh Sub53_Tw_3_derep_sorted.bam rg_out/Sub53_Tw_3.bam Sub53_Tw_3
bash addRG.sh UC_B11889_Env_2016_derep_sorted.bam rg_out/UC_B11889_Env_2016.bam UC_B11889_Env_2016
bash addRG.sh UC_B11891_Env_2016_derep_sorted.bam rg_out/UC_B11891_Env_2016.bam UC_B11891_Env_2016
swarm -f addrg.swarm   --job-name addrg -t 1 -g 98 --time 2:00:00 --logdir addRG_swarm_out  --gres=lscratch:400 

#index the bam files
cd /data/proctordm/cauris_genomes/alignments/rg_out/

bam_files=($(ls  *.bam ))
# get total subscripts in an array
total=${#bam_files[*]}

#loop over all the files in the input directories 
for (( i=0; i<=$(( $total -1 )); i++ ))
do
        samtools index ${bam_files[i]}
done

###### get the bam list
ls -d /data/proctordm/cauris_genomes/alignments/rg_out/*bam > bam.list



#### index the rerference genome
samtools faidx  Caur_3166.Fg_acuzp.v2.Nano.pilon.fasta
bwa index -p Caur.B11891_adqia.v2.Nano.pilon Caur.B11891_adqia.v2.Nano.pilon.fasta
java -Xmx96g -jar $PICARDJARPATH/picard.jar CreateSequenceDictionary -R Caur.B11891_adqia.v2.Nano.pilon.fasta



###################################################### call variants with freebayes parallel
#we will specify a quality score minimum of 30 since the default i think is 20 and a prior analysis revealed that this reduces number of variants significantly
###################################################### 
#!/usr/bin/bash
module load freebayes
module load bcftools
freebayes-parallel fb_regions.txt \
  36 \
  -f Caur.B11891_adqia.v2.Nano.pilon.fasta -p 1 \
  --min-base-quality 30 \
  -L bam.list \
  > freebayes_adqia_141samples_2022-09-29.vcf


bgzip freebayes_adqia_141samples_2022-09-29.vcf
tabix freebayes_adqia_141samples_2022-09-29.vcf.gz
bcftools view --min-ac 1:nref -Ov freebayes_adqia_141samples_2022-09-29.vcf.gz\
  > freebayes_adqia_141samples_2022-09-29_unsorted_unmasked.vcf


###################################################### call variants with bcftoools
#mpileup and then bcftools call
###################################################### 
#sbatch --mem=128g --cpus-per-task 32 --time 3-0 bcftools.sh
#!/usr/bin/bash
module load bcftools
#[+] Loading samtools 1.13 
bcftools mpileup -Ov -f /data/proctordm/cauris_genomes/reference/Caur.B11891_adqia.v2.Nano.pilon.fasta -b bam.list -B  -Q 30 | bcftools call -mv --ploidy 1 -o bcftools_adqia_141samples_2022-09-29.vcf


###################################################### call variants with gatk
#create a gatk haplotype caller on individual samples
###################################################### 

s
bam_names=($(ls  *.bam))
gatk_vcf_out=($(echo "${bam_names[*]}" | sed 's/.bam/.gatk_HC.g.vcf/g'))
total=${#bam_names[*]}
for (( i=0; i<=$(( $total -1 )); i++ ))
do	
echo "gatk --java-options '-Xmx8g' HaplotypeCaller -R /data/proctordm/cauris_genomes/reference/Caur.B11891_adqia.v2.Nano.pilon.fasta  -I ${bam_names[i]} -O ../gatk_out/${gatk_vcf_out[i]} -ploidy 1 -ERC GVCF"
done >gatk.swarm 


swarm -f gatk.swarm --job-name gatk -t 16 -g 96  --time 24:00:00 --module GATK --logdir ../gatk_out/gatk_swarm_out

#############################################joint genotyping
#nice reference: https://www.programmersought.com/article/10963095137/

#create a file called: intervals.list containing the following

########################group the vcsf so that we can perform joint genotyping

#create the intervals list
cat Caur.B11891_adqia.v2.Nano.pilon.fasta | grep '>' | grep -Eo '[^>]*' | grep -Eo '^[^ ]*' > intervals.list

#The Genome Analysis Toolkit (GATK) v4.2.4.1
#!/usr/bin/bash
module load GATKs
gatk GenomicsDBImport \
-V FL_B12847_urine_2019.gatk_HC.g.vcf \
-V FL_B12847_urine_2021.gatk_HC.g.vcf \
-V FL_B16329_blood_2020.gatk_HC.g.vcf \
-V IL_B11842_urine_2019.gatk_HC.g.vcf \
-V IL_B11843_blood_2019.gatk_HC.g.vcf \
-V IL_B11880_nares_2019.gatk_HC.g.vcf \
-V IL_B11882_nares_2019.gatk_HC.g.vcf \
-V IL_B11883_Rectal_2019.gatk_HC.g.vcf \
-V IL_B11884_Vaginal_2019.gatk_HC.g.vcf \
-V IL_B11885_Axilla_2019.gatk_HC.g.vcf \
-V IL_B11886_Groin_2019.gatk_HC.g.vcf \
-V IL_B12028_blood_2019.gatk_HC.g.vcf \
-V IL_B12029_blood_2019.gatk_HC.g.vcf \
-V IL_B12030_blood_2019.gatk_HC.g.vcf \
-V IL_B12031_blood_2019.gatk_HC.g.vcf \
-V IL_B12032_blood_2019.gatk_HC.g.vcf \
-V IL_B12033_blood_2019.gatk_HC.g.vcf \
-V IL_B12034_blood_2019.gatk_HC.g.vcf \
-V IL_B12035_blood_2019.gatk_HC.g.vcf \
-V IL_B12036_blood_2019.gatk_HC.g.vcf \
-V IL_B12046_AxillaGroin_2019.gatk_HC.g.vcf \
-V IL_B12077_feces_2019.gatk_HC.g.vcf \
-V IL_B12189_nares_2019.gatk_HC.g.vcf \
-V IL_B12190_Axilla_2019.gatk_HC.g.vcf \
-V IL_B12229_kidney_2019.gatk_HC.g.vcf \
-V IL_B12378_AxillaGroin_2019.gatk_HC.g.vcf \
-V IL_B12380_AxillaGroin_2019.gatk_HC.g.vcf \
-V IL_B12388_AxillaGroin_2019.gatk_HC.g.vcf \
-V IL_B12406_urine_2019.gatk_HC.g.vcf \
-V IL_B12510_AxillaGroin_2019.gatk_HC.g.vcf \
-V IL_B12579_Axilla_2019.gatk_HC.g.vcf \
-V IL_B12937_Axilla_2019.gatk_HC.g.vcf \
-V IL_B12938_Axilla_2019.gatk_HC.g.vcf \
-V IL_CA01_Pleuralfluid_2021.gatk_HC.g.vcf \
-V IL_CA03_urine_2021.gatk_HC.g.vcf \
-V IL_CA04_respiratory_2021.gatk_HC.g.vcf \
-V IL_CA05_blood_2021.gatk_HC.g.vcf \
-V IL_CA06_blood_2021.gatk_HC.g.vcf \
-V IL_CA07_wound_2021.gatk_HC.g.vcf \
-V IL_CA08_Patientroom_2021.gatk_HC.g.vcf \
-V IL_CA09_Patientroom_2021.gatk_HC.g.vcf \
-V IL_CA10_Patientroom_2021.gatk_HC.g.vcf \
-V IL_CA11_urine_2021.gatk_HC.g.vcf \
-V IL_CA13_wound_2021.gatk_HC.g.vcf \
-V IL_CA14_urine_2021.gatk_HC.g.vcf \
-V IL_CA15_skin_2021.gatk_HC.g.vcf \
-V IL_CA16_Nose_2021.gatk_HC.g.vcf \
-V IL_CA17_Axilla_2021.gatk_HC.g.vcf \
-V IL_CA18_Groin_2021.gatk_HC.g.vcf \
-V IL_CA19_skin_2021.gatk_HC.g.vcf \
-V IL_CA20_skin_2021.gatk_HC.g.vcf \
-V IL_CA21_respiratory_2021.gatk_HC.g.vcf \
-V IL_CA22_blood_2021.gatk_HC.g.vcf \
-V IL_CA23_urine_2021.gatk_HC.g.vcf \
-V IL_CA24_urine_2021.gatk_HC.g.vcf \
-V IL_CA25_Synovialfluid_2021.gatk_HC.g.vcf \
-V K_15160_Urine_2018.gatk_HC.g.vcf \
-V K_B12189_Nares_2016.gatk_HC.g.vcf \
-V K_B12510_AxGroin_2017.gatk_HC.g.vcf \
-V K_B14707_AxGroin_2018.gatk_HC.g.vcf \
-V K_B14714_AxGroin_2018.gatk_HC.g.vcf \
-V K_B15161_Blood_2018.gatk_HC.g.vcf \
-V K_B16276_Blood_2018.gatk_HC.g.vcf \
-V K_B16278_Blood_2018.gatk_HC.g.vcf \
-V MA_B12493_Bronchus_2019.gatk_HC.g.vcf \
-V Sub14_An_1.gatk_HC.g.vcf \
-V Sub14_An_2.gatk_HC.g.vcf \
-V Sub14_Fg_1.gatk_HC.g.vcf \
-V Sub14_Fg_2.gatk_HC.g.vcf \
-V Sub14_Ic_3.gatk_HC.g.vcf \
-V Sub14_N_1.gatk_HC.g.vcf \
-V Sub14_N_2.gatk_HC.g.vcf \
-V Sub14_N_3.gatk_HC.g.vcf \
-V Sub14_Tw_1.gatk_HC.g.vcf \
-V Sub15_An_1.gatk_HC.g.vcf \
-V Sub15_Fg_2.gatk_HC.g.vcf \
-V Sub15_Fg_3.gatk_HC.g.vcf \
-V Sub15_Ic_1.gatk_HC.g.vcf \
-V Sub15_Tw_2.gatk_HC.g.vcf \
-V Sub15_Tw_3.gatk_HC.g.vcf \
-V Sub2_Ax_1.gatk_HC.g.vcf \
-V Sub2_Fg_1.gatk_HC.g.vcf \
-V Sub2_Fg_2.gatk_HC.g.vcf \
-V Sub2_Fg_3.gatk_HC.g.vcf \
-V Sub2_Ic_2.gatk_HC.g.vcf \
-V Sub2_N_2.gatk_HC.g.vcf \
-V Sub2_N_3.gatk_HC.g.vcf \
-V Sub23_An_1.gatk_HC.g.vcf \
-V Sub23_Ea_2.gatk_HC.g.vcf \
-V Sub23_Ea_3.gatk_HC.g.vcf \
-V Sub23_Fg_1.gatk_HC.g.vcf \
-V Sub23_Ic_3.gatk_HC.g.vcf \
-V Sub23_N_1.gatk_HC.g.vcf \
-V Sub23_N_2.gatk_HC.g.vcf \
-V Sub23_N_3.gatk_HC.g.vcf \
-V Sub28_Fg_1.gatk_HC.g.vcf \
-V Sub28_Fg_3.gatk_HC.g.vcf \
-V Sub28_N_1.gatk_HC.g.vcf \
-V Sub28_N_2.gatk_HC.g.vcf \
-V Sub28_N_3.gatk_HC.g.vcf \
-V Sub35_An_3.gatk_HC.g.vcf \
-V Sub35_Ic_2.gatk_HC.g.vcf \
-V Sub35_Ic_3.gatk_HC.g.vcf \
-V Sub35_N_2.gatk_HC.g.vcf \
-V Sub35_N_3.gatk_HC.g.vcf \
-V Sub4_Ic_1.gatk_HC.g.vcf \
-V Sub4_Ic_3.gatk_HC.g.vcf \
-V Sub4_N_1.gatk_HC.g.vcf \
-V Sub4_N_2.gatk_HC.g.vcf \
-V Sub4_N_3.gatk_HC.g.vcf \
-V Sub4_Tw_2.gatk_HC.g.vcf \
-V Sub46_An_1.gatk_HC.g.vcf \
-V Sub46_Fg_1.gatk_HC.g.vcf \
-V Sub46_Fg_2.gatk_HC.g.vcf \
-V Sub46_Fg_3.gatk_HC.g.vcf \
-V Sub46_N_1.gatk_HC.g.vcf \
-V Sub46_N_2.gatk_HC.g.vcf \
-V Sub46_N_3.gatk_HC.g.vcf \
-V Sub46_Ne_2.gatk_HC.g.vcf \
-V Sub46_Ne_3.gatk_HC.g.vcf \
-V Sub48_An_1.gatk_HC.g.vcf \
-V Sub48_An_3.gatk_HC.g.vcf \
-V Sub48_Fg_1.gatk_HC.g.vcf \
-V Sub48_Fg_3.gatk_HC.g.vcf \
-V Sub48_Ic_1.gatk_HC.g.vcf \
-V Sub48_N_1.gatk_HC.g.vcf \
-V Sub48_N_3.gatk_HC.g.vcf \
-V Sub48_Tw_1.gatk_HC.g.vcf \
-V Sub5_An_1.gatk_HC.g.vcf \
-V Sub5_Ea_2.gatk_HC.g.vcf \
-V Sub5_N_1.gatk_HC.g.vcf \
-V Sub5_N_3.gatk_HC.g.vcf \
-V Sub53_Fg_1.gatk_HC.g.vcf \
-V Sub53_Fg_3.gatk_HC.g.vcf \
-V Sub53_N_1.gatk_HC.g.vcf \
-V Sub53_N_2.gatk_HC.g.vcf \
-V Sub53_N_3.gatk_HC.g.vcf \
-V Sub53_Tw_2.gatk_HC.g.vcf \
-V Sub53_Tw_3.gatk_HC.g.vcf \
-V UC_B11889_Env_2016.gatk_HC.g.vcf \
-V UC_B11891_Env_2016.gatk_HC.g.vcf \
-L /data/proctordm/cauris_genomes/reference/intervals.list \
    --genomicsdb-workspace-path my_database 

#cleanup
mkdir  mkdir haplotype_caller_out
	mv *gatk_HC*  mkdir haplotype_caller_out
	
############ GenotypeGVCFs
gatk --java-options "-Xmx4g" GenotypeGVCFs \
   -R /data/proctordm/cauris_genomes/reference/Caur.B11891_adqia.v2.Nano.pilon.fasta\
   -V gendb://my_database \
   -O gatk141.GenotypeGVCFs.vcf.gz
 
 
#################################### 
#filter and mask
##################################### 


#mask the vcf files based on repeat masker
bedtools intersect -v -a bcftools_adqia_141samples_2022-09-29.vcf  -b /data/proctordm/cauris_genomes/0.genome_prep_out/mask/repeat_masker.bed -wa -header > bcftools_masked.vcf 
bedtools intersect -v -a freebayes_adqia_141samples_2022-09-29_unsorted_unmasked.vcf  -b /data/proctordm/cauris_genomes/0.genome_prep_out/mask/repeat_masker.bed -wa -header > freebayes_masked.vcf 
bedtools intersect -v -a gatk141.GenotypeGVCFs.vcf.gz  -b /data/proctordm/cauris_genomes/0.genome_prep_out/mask/repeat_masker.bed -wa -header > gatk_masked.vcf 

mv * ../6.masked



#Add all possible tags to vcf file
for file in *.vcf; do
  bcftools +fill-tags $file -- -t ALL > tmp
  mv tmp $file
done

#get the number of samples
bcftools query -l bcftools_masked.vcf| sort > samples.txt

mkdir -p unfiltered
for file in *_masked.vcf; do
  cp $file $(pwd)/unfiltered/
done

# filter the SNPs based on three filters:
# min. 5 reads per sample
# minimum quality score 1
# only biallelic SNPs

num_samples=$(wc -l < samples.txt)
depth_filter=$((num_samples * 5))

bcftools view -v snps bcftools_masked.vcf  | bcftools filter -i 'MIN(INFO/DP)>'${depth_filter} | bcftools filter -i 'QUAl>30' | bcftools view -m2 -M2 -v snps  > bcftools_filtered_masked.vcf
bcftools view -v snps freebayes_masked.vcf| bcftools filter -i 'MIN(INFO/DP)>'${depth_filter} | bcftools filter -i 'QUAl>30' | bcftools view -m2 -M2 -v snps  > freebayes_filtered_masked.vcf
bcftools view -v snps 	gatk_masked.vcf | bcftools filter -i 'MIN(INFO/DP)>'${depth_filter} | bcftools filter -i 'QUAl>30' | bcftools view -m2 -M2 -v snps  > gatk_filtered_masked.vcf







###########################################
#combine calllers
##########################################

# get the sample names and put them into a file so we can sort freebayes
# and gatk in the same order before combining
bcftools query -l bcftools_filtered_masked.vcf  | sort > samples.txt

# sort the samples of all three outputs
bcftools view -S samples.txt bcftools_unmasked.vcf \
  > bcf_tmp
mv bcf_tmp bcftools_unmasked.vcf

# if we're working with the vcf files where the output
bcftools view -S samples.txt freebayes_filtered_masked.vcf > freebayes_filtered_masked_sorted.vcf 
bcftools view -S samples.txt gatk_filtered_masked.vcf > gatk_filtered_masked_sorted.vcf 

module load GATK

# # normalize variants
# # https://genome.sph.umich.edu/wiki/Variant_Normalization
echo "[$(date)] ${fname}: Normalizing variants." >> $log_file
for vcf in bcftools_filtered_masked.vcf  freebayes_filtered_masked_sorted.vcf  gatk_filtered_masked_sorted.vcf; do
    bcftools norm -f /data/proctordm/cauris_genomes/reference/Caur.B11891_adqia.v2.Nano.pilon.fasta $vcf -Ov > tmp
    mv tmp $vcf
done
	#Lines   total/split/realigned/skipped:	909/0/0/0
	#Lines   total/split/realigned/skipped:	823/0/20/0
	#Lines   total/split/realigned/skipped:	951/0/0/0


# we need to zip and index the files before passing them to bcftools
 

bgzip bcftools_filtered_masked.vcf  
bgzip freebayes_filtered_masked_sorted.vcf 
bgzip gatk_filtered_masked_sorted.vcf
tabix bcftools_filtered_masked.vcf.gz
tabix freebayes_filtered_masked_sorted.vcf.gz
tabix gatk_filtered_masked_sorted.vcf.gz

# take the 3-way intersection of all variant callers
bcftools isec -c all -n=3 -Ov -p 3way\
  bcftools_filtered_masked.vcf.gz \
  freebayes_filtered_masked_sorted.vcf.gz \
  gatk_filtered_masked_sorted.vcf.gz

# bcftools isec produces a lot of extra files that we
# don't actually want, so we delete them
rm 0000.vcf
mv 0001.vcf 3_way_masked.vcf
rm 0002.vcf
rm README.txt
rm sites.txt


#################
#select variants for the one sample that should be a control
bcftools view --samples-file samples2keep.txt --min-ac=1 --no-update 3_way_masked.vcf.gz > adqia.3way.bcftools.vcf
#control the adqia vcf to a bed file format
sed -e 's/chr//' adqia.3way.bcftools.vcf| awk '{OFS="\t"; if (!/^#/){print $1,$2-1,$2,$4"/"$5,"+"}}'> adqia.bed 
#remove the SNPs that were present in the negative control
bedtools intersect -v -a 3_way_masked.vcf.gz  -b adqia.bed  -wa -header > 3_way_masked_adqia_masked.vcf


###################################################### 
# convert to a phylip object and make a tree
###################################################### 
#two directories her
#1: contains the shared variant consensus files: /data/proctordm/comp_population/population71/consensus
#2: contains the variant caller specific files: /data/proctordm/comp_population/population71/varcaller
#########now let's convert the vcf to a phylip file so we can generate phylogenetic trees from the freebayes data
mkdir vcf2phylip
cd vcf2phylip
cp /data/Segrelab/loading_dock/vcf2phylip/* .

### this script expects file endings to be either .gz or .rf so we need to do some cleanup
cd /data/proctordm/cauris_genomes/8.3wayconsensus
#check file type
file *vcf

bgzip 3_way_masked_adqia_masked.vcf
tabix 3_way_masked_adqia_masked.vcf.gz

#convert to phylip
python /data/Segrelab/loading_dock/vcf2phylip/vcf2phylip.py --input 3_way_masked_adqia_masked.vcf.gz

#########now let's generate trees
phy1=3_way_masked_adqia_masked.min4.phy

module load iqtree
#[+] Loading iqtree  2.1.2

python /data/Segrelab/loading_dock/vcf2phylip/vcf2phylip.py --input 3_way_masked_adqia_masked_snpeff_OutbreakOnly.vcf 

bbtools phylip2fasta in=3_way_masked_adqia_masked_snpeff_OutbreakOnly.min4.phy out=3_way_masked_adqia_masked_snpeff_OutbreakOnly.min4.fasta
ADQIA_annot.gtf 

split -d -l 10 ADQIA_annot.gtf chr3_chunks
####### get trees for the variant caller specific files
iqtree2 -s $phy1  -T AUTO -B 10000
	NOTE: FL_B12847_urine_2021 is identical to FL_B12847_urine_2019 but kept for subsequent analysis
	NOTE: IL_B12229_kidney_2019 is identical to IL_B11842_urine_2019 but kept for subsequent analysis
	NOTE: IL_B11882_nares_2019 is identical to IL_B11880_nares_2019 but kept for subsequent analysis
	NOTE: IL_B12388_AxillaGroin_2019 is identical to IL_B12046_AxillaGroin_2019 but kept for subsequent analysis
	NOTE: IL_B12190_Axilla_2019 is identical to IL_B12189_nares_2019 but kept for subsequent analysis
	NOTE: K_B12510_AxGroin_2017 is identical to IL_B12510_AxillaGroin_2019 but kept for subsequent analysis
	NOTE: Sub15_An_1 is identical to Sub14_An_1 but kept for subsequent analysis
	NOTE: Sub15_Ic_1 is identical to Sub15_Fg_2 but kept for subsequent analysis
	NOTE: Sub23_Ea_3 is identical to Sub23_Ea_2 but kept for subsequent analysis
	NOTE: Sub2_Fg_2 is identical to Sub2_Ax_1 but kept for subsequent analysis
	NOTE: Sub48_Fg_3 is identical to Sub48_An_3 but kept for subsequent analysis
	NOTE: Sub4_N_3 is identical to Sub4_N_1 but kept for subsequent analysis
	NOTE: Sub53_N_2 is identical to Sub53_N_1 but kept for subsequent analysis
	NOTE: 13 identical sequences (see below) will be ignored for subsequent analysis
	NOTE: IL_B11885_Axilla_2019 (identical to IL_B11880_nares_2019) is ignored but added at the end
	NOTE: IL_B12031_blood_2019 (identical to IL_B11880_nares_2019) is ignored but added at the end
	NOTE: IL_B12032_blood_2019 (identical to IL_B11880_nares_2019) is ignored but added at the end
	NOTE: IL_B12033_blood_2019 (identical to IL_B11880_nares_2019) is ignored but added at the end
	NOTE: IL_B12034_blood_2019 (identical to IL_B11880_nares_2019) is ignored but added at the end
	NOTE: IL_B12035_blood_2019 (identical to IL_B11880_nares_2019) is ignored but added at the end
	NOTE: IL_B12036_blood_2019 (identical to IL_B11880_nares_2019) is ignored but added at the end
	NOTE: UC_B11889_Env_2016 (identical to IL_B11880_nares_2019) is ignored but added at the end
	NOTE: IL_CA24_urine_2021 (identical to IL_B12046_AxillaGroin_2019) is ignored but added at the end
	NOTE: IL_B12378_AxillaGroin_2019 (identical to IL_B12189_nares_2019) is ignored but added at the end
	NOTE: K_B12189_Nares_2016 (identical to IL_B12189_nares_2019) is ignored but added at the end
	NOTE: Sub2_Fg_3 (identical to Sub2_Ax_1) is ignored but added at the end
	NOTE: Sub2_N_2 (identical to Sub2_Ax_1) is ignored but added at the end



###################################################### 
# Get functional impacts using snpeff
###################################################### 

#How to create a snpEff database using a gff3 and genomic DNA fasta file... (note, the chromosome names must match in the 2 files)

DBNAME='Cauris_CladeIV_ADQIA'
GFF3='/data/proctordm/cauris_genomes/reference/Caur.B11891_adqia.v2.Nano.pilon.ff.gff'
FASTA='/data/proctordm/cauris_genomes/reference/Caur.B11891_adqia.v2.Nano.pilon.fasta'

 #Go into the snpEff directory and create a directory for your files
cd /home/proctordm/bin/snpEff
mkdir data
mkdir data/$DBNAME

 #Copy the files into snpEff's directory structure
cp $GFF3 data/$DBNAME/genes.gff
cp $FASTA data/$DBNAME/sequences.fa

 #Edit snpEff.config and insert your specific database information:
echo "$DBNAME.genome : $DBNAME" >> snpEff.config

 #Build the database
 java -jar snpEff.jar build -gff3 -v $DBNAME
 
 

#annotate the variants 
cd /home/proctordm/bin/snpEff
VCF='/data/proctordm/cauris_genomes/9.input/3_way_masked_adqia_masked.vcf'
java -Xmx8g -jar  /home/proctordm/bin/snpEff/snpEff.jar -c snpEff.config $DBNAME $VCF   > $VCF_sneff_out1
	mkdir consensus3way
		mv snpEff_summary.html consensus3way/
		mv snpEff_genes.txt consensus3way/
		mv $VCF_sneff_out1 consensus3way


java -Xmx8g -jar /home/proctordm/bin/snpEff/snpEff.jar -c snpEff.config $DBNAME gatk.sharedPositions2way_71genomes_bcftools_freebayes.SNP.population.vcf.gz > gatk.sharedPositions2way_snpeff.vcf
	mkdir consensus2way
		mv snpEff_summary.html consensus2way/
		mv snpEff_genes.txt consensus2way/
		mv gatk.sharedPositions2way_snpeff.vcf consensus2way


 mkdir snpeff_out
 	mv consensus3way snpeff_out
  	mv consensus2way snpeff_out
  	
  	
 

#####################################################
#use esearrch to get the gene names and fasta sequences for eachh of the significant hits
#####################################################
module load edirect
#[+] Loading edirect  14.5 
cd /data/proctordm/comp_population/population71/consensus/snpeff_out/consensus2way
awk '{print $2}' snpEff_genes.txt > mygenes.2way


#use nano to get rid of the first two linees
readarray myArray <  mygenes.2way
total=${#myArray[*]}
for (( i=0; i<=$(( $total -1 )); i++ ))
do
esearch -db protein -query ${myArray[i]} | efetch -format fasta > ${myArray[i]}
done

head -n1 -q CJ* > fasta.headers.2way.txt
  	
  	


#############################
# summarize functional predictions
#############################
cd /data/proctordm/comp_population/population71/consensus/snpeff_out/consensus2way

#use snpsift to get a one by one linee summarry of the impacts
#java -Xmx8g -jar /home/proctordm/bin/snpEff/SnpSift.jar extractFields gatk.sharedPositions2way_snpeff.vcf CHROM POS REF ALT "ANN[*].EFFECT" "ANN[*].GENEID" > gatk.sharedPositions2way_snpeff_snpsift.txt
#java -Xmx8g -jar /home/proctordm/bin/snpEff/SnpSift.jar extractFields -s "," -e "." gatk.sharedPositions2way_snpeff.vcf  CHROM POS REF ALT  "EFF[*].EFFECT" "EFF[*].AA" "ANN[*].EFFECT" "ANN[*].GENEID" >  gatk.sharedPositions2way_snpeff_snpsift_onebyonee.txt
cat gatk.sharedPositions2way_snpeff.vcf | perl /home/proctordm/bin/snpEff/scripts/vcfEffOnePerLine.pl | java -Xmx8g -jar /home/proctordm/bin/snpEff/SnpSift.jar extractFields  - CHROM POS REF ALT AF "EFF[*].EFFECT" "EFF[*].AA" "ANN[*].EFFECT" "ANN[*].GENEID" "ANN[*].BIOTYPE" "ANN[*].CODON" "ANN[*].IMPACT" >  gatk.sharedPositions2way_snpeff_snpsift_onebyone_AF.txt


#us vcf2tabl to get a nice summary of the variants
java -jar /home/proctordm/jvarkit/dist/vcf2table.jar vcfpredictions_results.vcf --color --format html > out.html
java -jar /home/proctordm/jvarkit/dist/vcf2table.jar gatk.sharedPositions2way_snpeff.vcf --color --format text > snepeff_vcf2table.out.txt

#use vcf2table to get predictioons also 
java -Xmx4g -jar $PICARDJARPATH/picard.jar  CreateSequenceDictionary R=/data/proctordm/comp_population/population71/consensus/data/Cauris_CladeIV_ACUZR/sequences.fa O=/data/proctordm/comp_population/population71/consensus/data/Cauris_CladeIV_ACUZR/sequences.dict
java -jar /home/proctordm/jvarkit/dist/vcf2table.jar gatk.sharedPositions2way_snpeff.vcf --color --format html > out.html
java -Xmx4g -jar /home/proctordm/jvarkit/dist/vcfpredictions.jar  gatk.sharedPositions2way_snpeff.vcf -k /data/proctordm/comp_population/population71/consensus/data/Cauris_CladeIV_ACUZR/genes.gtf -R /data/proctordm/comp_population/population71/consensus/data/Cauris_CladeIV_ACUZR/sequences.fa > vcfpredictions_results.vcf

java -jar /home/proctordm/jvarkit/dist/vcf2table.jar vcfpredictions_results.vcf  --format text > vcfpredictions_results_vcf2table_out.txt

#https://www.biostars.org/p/298361/
./gradlew bioalcidaejdk
java -jar /home/proctordm/jvarkit/dist/bioalcidaejdk.jar -e 'final List<GenotypeType> all_types = Arrays.asList(GenotypeType.values());println("CHROM\tPOS\tREF\tALT\tAF\t"+all_types.stream().map(T->"count."+T.name()+"\tsample."+T.name()).collect(Collectors.joining("\t")));stream().forEach(V->println( V.getContig()+"\t"+ V.getStart()+"\t"+ V.getReference().getDisplayString()+"\t"+ V.getAlternateAlleles().stream().map(A->A.getDisplayString()).collect(Collectors.joining(","))+"\t"+ all_types.stream().map(T-> String.valueOf( V.getGenotypes().stream(). filter(G->G.getType().equals(T)). count()) + "\t"+ V.getGenotypes().stream(). filter(G->G.getType().equals(T)). map(G->G.getSampleName()). collect(Collectors.joining(",")) ). collect(Collectors.joining("\t"))) );'   gatk.sharedPositions2way_71genomes_gatk_freebayes.SNP.population_snpeff.vcf >  gatk.sharedPositions2way_metagenome.SNP.population_snpeff_out_samplesummarry.vcf 


########################
#get snps private to each vcf
####c#################
cd /data/proctordm/comp_population/subjectwise

######## bgzip files
vcf=($(ls  *.vcf))
# get total subscripts in an array
	total=${#vcf[*]}
	#loop over all the files in the input directories 
	for (( i=0; i<=$(( $total -1 )); i++ ))
	do
		bgzip ${vcf[i]}
	done

######## tabix files
vcf=($(ls  *.gz))
# get total subscripts in an array
	total=${#vcf[*]}
	#loop over all the files in the input directories 
	for (( i=0; i<=$(( $total -1 )); i++ ))
	do
		tabix --preset vcf ${vcf[i]}
	done

######## intersection
load bedtools
#[+] Loading bedtools  2.30.0

subject23='gatk71.subject23.GenotypeGVCFs.vcf.gz'  
subject28='gatk71.subject28.GenotypeGVCFs.vcf.gz'
subject46='gatk71.subject46.GenotypeGVCFs.vcf.gz'
subject53='gatk71.subject53.GenotypeGVCFs.vcf.gz'

bcftools view -v snps $subject23 | grep -v "^#" | wc -l #458
bcftools view -v snps $subject28 | grep -v "^#" | wc -l #414
bcftools view -v snps $subject46 | grep -v "^#" | wc -l #512
bcftools view -v snps $subject53 | grep -v "^#" | wc -l #460



bcftools isec -p allway -n=1 -c all $subject23 $subject28 $subject35 $subject46 $subject53
bcftools isec -p sub23vsub28 -n=1  $subject23 $subject28 
bcftools isec -p sub23vsub35 -n=1 $subject23 $subject35
bcftools isec -p sub23vsub46 -n=1 $subject23 $subject46
bcftools isec -p sub23vsub53 -n=1 $subject23 $subject53

bcftools isec -p sub28vsub35 -n=1 $subject28 $subject35
bcftools isec -p sub28vsub46 -n=1 $subject28 $subject46
bcftools isec -p sub28vsub53 -n=1 $subject28 $subject53

bcftools isec -p sub35vsub46 -n=1 $subject46 $subject35


subject23=gatk.sharedPositions2way_Subject23.SNP.population.vcf.gz
subject28=gatk.sharedPositions2way_Subject28.SNP.population.vcf.gz
subject35=gatk.sharedPositions2way_Subject35.SNP.population.vcf.gz
subject46=gatk.sharedPositions2way_Subject46.SNP.population.vcf.gz
subject53=gatk.sharedPositions2way_Subject53.SNP.population.vcf.gz
tabix --preset vcf $subject23
tabix --preset vcf $subject28
tabix --preset vcf $subject35
tabix --preset vcf $subject46
tabix --preset vcf $subject53
bcftools isec -p allway -n=1 -c all $subject23 $subject28 $subject35 $subject46 $subject53
bcftools isec -p sub23vsub28 -n=1  $subject23 $subject28 
bcftools isec -p sub23vsub35 -n=1 $subject23 $subject35
bcftools isec -p sub23vsub46 -n=1 $subject23 $subject46
bcftools isec -p sub23vsub53 -n=1 $subject23 $subject53

bcftools isec -p sub28vsub35 -n=1 $subject28 $subject35
bcftools isec -p sub28vsub46 -n=1 $subject28 $subject46
bcftools isec -p sub28vsub53 -n=1 $subject28 $subject53

bcftools isec -p sub35vsub46 -n=1 $subject46 $subject35


###################################################### 
# Get coverrage
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


mkdir bbtools_out
	mv *.bbtools.coverage.txt bbtools_out

mkdir swarm_out
	mv *34281740* swarm_out/


