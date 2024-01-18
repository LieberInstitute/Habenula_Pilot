snRNA-seq FASTQs
================

From scripts at https://github.com/LieberInstitute/Habenula_Pilot/tree/master/code/07_cellranger such as https://github.com/LieberInstitute/Habenula_Pilot/blob/master/code/07_cellranger/jobSubmit_21-3-22_enelson_Br1204-Hb.sh, we can locate the `--fastqs` and `--sample` tags to find the original file locations when these scripts were run.

```
/dcl01/ajaffe/data/lab/singleCell/10x_pilot/FASTQ/Feb2021/Br1204_Hb/Br1204-Hb
/dcl01/ajaffe/data/lab/singleCell/10x_pilot/FASTQ/Feb2021/Br5558_Hb/Br5558-Hb
/dcl01/ajaffe/data/lab/singleCell/10x_pilot/FASTQ/Br1469_Hb/Br1469_Hb
/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/data_from_linda/2021-05-21/ASpa050421/3c-k
/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/data_from_linda/2021-05-21/ASpa050421/1c-k
/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/data_from_linda/2021-05-21/ASpa050421/4c-k
/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/data_from_linda/2021-05-21/ASpa050421/2c-k
```

`/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/data_from_linda/2021-05-21/ASpa050421` was moved to `/dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2021-05-21_ASpa050421/` and other files were moved to `/dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2021-03-15_KMay022421/` and `/dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2021-10-15_KMay110220/`. (Internal note: check this [Slack DM thread](https://jhu-genomics.slack.com/archives/C06EDFXF03C/p1705533086510189) for more details.)

Thus the files are now located at:

```
/dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2021-03-15_KMay022421/Br1204_Hb* .
/dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2021-03-15_KMay022421/Br5558_Hb*
/dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2021-10-15_KMay110220/Br1469_Hb*
/dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2021-05-21_ASpa050421/3c_k*
/dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2021-05-21_ASpa050421/1c_k*
/dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2021-05-21_ASpa050421/4c_k*
/dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2021-05-21_ASpa050421/2c_k*
```

as we can see by listing their contents:

```bash
ls -lh /dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2021-03-15_KMay022421/Br1204_Hb*
# /dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2021-03-15_KMay022421/Br1204_Hb_L001_ds.450f8b4b14ce45d99a81d702fad67d74:
# total 6.4G
# -rwxrwx---+ 1 lcollado lieber_lmh_storage 2.0G Jun 17  2021 Br1204-Hb_S3_L001_R1_001.fastq.gz
# -rwxrwx---+ 1 lcollado lieber_lmh_storage 4.4G Jun 17  2021 Br1204-Hb_S3_L001_R2_001.fastq.gz
# 
# /dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2021-03-15_KMay022421/Br1204_Hb_L002_ds.3ebb12ef36b041faba06599157a53d17:
# total 5.9G
# -rwxrwx---+ 1 lcollado lieber_lmh_storage 1.9G Jun 17  2021 Br1204-Hb_S3_L002_R1_001.fastq.gz
# -rwxrwx---+ 1 lcollado lieber_lmh_storage 4.1G Jun 17  2021 Br1204-Hb_S3_L002_R2_001.fastq.gz


ls -lh /dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2021-03-15_KMay022421/Br5558_Hb*
# /dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2021-03-15_KMay022421/Br5558_Hb_L001_ds.68bb30be674d4e3da8abd44466ecf985:
# total 7.0G
# -rwxrwx---+ 1 lcollado lieber_lmh_storage 2.2G Jun 17  2021 Br5558-Hb_S2_L001_R1_001.fastq.gz
# -rwxrwx---+ 1 lcollado lieber_lmh_storage 4.9G Jun 17  2021 Br5558-Hb_S2_L001_R2_001.fastq.gz
# 
# /dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2021-03-15_KMay022421/Br5558_Hb_L002_ds.c7f8800daced4fddaee43c22fe006a4c:
# total 6.7G
# -rwxrwx---+ 1 lcollado lieber_lmh_storage 2.1G Jun 17  2021 Br5558-Hb_S2_L002_R1_001.fastq.gz
# -rwxrwx---+ 1 lcollado lieber_lmh_storage 4.7G Jun 17  2021 Br5558-Hb_S2_L002_R2_001.fastq.gz

ls -lh /dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2021-10-15_KMay110220/Br1469_Hb*
# total 14G
# -rwxrwx---+ 1 lcollado lieber_lmh_storage 2.1G Nov 16  2020 Br1469_Hb_S13_L001_R1_001.fastq.gz
# -rwxrwx---+ 1 lcollado lieber_lmh_storage 4.6G Nov 16  2020 Br1469_Hb_S13_L001_R2_001.fastq.gz
# -rwxrwx---+ 1 lcollado lieber_lmh_storage 2.1G Nov 16  2020 Br1469_Hb_S13_L002_R1_001.fastq.gz
# -rwxrwx---+ 1 lcollado lieber_lmh_storage 4.7G Nov 16  2020 Br1469_Hb_S13_L002_R2_001.fastq.gz

ls -lh /dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2021-05-21_ASpa050421/3c_k*
# /dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2021-05-21_ASpa050421/3c_k_L001_ds.19d0c09d1d5b4922add0828bf38422e2:
# total 5.2G
# -rwxrwx---+ 1 lcollado lieber_lmh_storage 1.6G Jun 17  2021 3c-k_S15_L001_R1_001.fastq.gz
# -rwxrwx---+ 1 lcollado lieber_lmh_storage 3.7G Jun 17  2021 3c-k_S15_L001_R2_001.fastq.gz
# 
# /dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2021-05-21_ASpa050421/3c_k_L002_ds.ee8551ea791a4794ba2553444b7e88fe:
# total 5.2G
# -rwxrwx---+ 1 lcollado lieber_lmh_storage 1.6G Jun 17  2021 3c-k_S15_L002_R1_001.fastq.gz
# -rwxrwx---+ 1 lcollado lieber_lmh_storage 3.7G Jun 17  2021 3c-k_S15_L002_R2_001.fastq.gz

ls -lh /dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2021-05-21_ASpa050421/1c_k*
# /dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2021-05-21_ASpa050421/1c_k_L001_ds.0a7b7047146a4e8580777073634c29ae:
# total 5.9G
# -rwxrwx---+ 1 lcollado lieber_lmh_storage 1.9G Jun 17  2021 1c-k_S13_L001_R1_001.fastq.gz
# -rwxrwx---+ 1 lcollado lieber_lmh_storage 4.0G Jun 17  2021 1c-k_S13_L001_R2_001.fastq.gz
# 
# /dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2021-05-21_ASpa050421/1c_k_L002_ds.ed2439f840eb47d08d139e21483725f6:
# total 5.9G
# -rwxrwx---+ 1 lcollado lieber_lmh_storage 1.9G Jun 17  2021 1c-k_S13_L002_R1_001.fastq.gz
# -rwxrwx---+ 1 lcollado lieber_lmh_storage 4.0G Jun 17  2021 1c-k_S13_L002_R2_001.fastq.gz

ls -lh /dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2021-05-21_ASpa050421/4c_k*
# /dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2021-05-21_ASpa050421/4c_k_L001_ds.1e5dca4b37b6461f9e6a2d460d532323:
# total 5.3G
# -rwxrwx---+ 1 lcollado lieber_lmh_storage 1.7G Jun 17  2021 4c-k_S16_L001_R1_001.fastq.gz
# -rwxrwx---+ 1 lcollado lieber_lmh_storage 3.6G Jun 17  2021 4c-k_S16_L001_R2_001.fastq.gz
# 
# /dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2021-05-21_ASpa050421/4c_k_L002_ds.d5822165012f40008c02444152516e33:
# total 5.3G
# -rwxrwx---+ 1 lcollado lieber_lmh_storage 1.7G Jun 17  2021 4c-k_S16_L002_R1_001.fastq.gz
# -rwxrwx---+ 1 lcollado lieber_lmh_storage 3.6G Jun 17  2021 4c-k_S16_L002_R2_001.fastq.gz

ls -lh /dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2021-05-21_ASpa050421/2c_k*
# /dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2021-05-21_ASpa050421/2c_k_L001_ds.01c025c8a9274af583a146419a7b973b:
# total 5.4G
# -rwxrwx---+ 1 lcollado lieber_lmh_storage 1.7G Jun 17  2021 2c-k_S14_L001_R1_001.fastq.gz
# -rwxrwx---+ 1 lcollado lieber_lmh_storage 3.8G Jun 17  2021 2c-k_S14_L001_R2_001.fastq.gz
# 
# /dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2021-05-21_ASpa050421/2c_k_L002_ds.c77782c03aef40769707e63d06135d8d:
# total 5.4G
# -rwxrwx---+ 1 lcollado lieber_lmh_storage 1.7G Jun 17  2021 2c-k_S14_L002_R1_001.fastq.gz
# -rwxrwx---+ 1 lcollado lieber_lmh_storage 3.8G Jun 17  2021 2c-k_S14_L002_R2_001.fastq.gz
```

We copied the files over with these commands:

```bash
pwd
# /dcs04/lieber/lcolladotor/pilotHb_LIBD001/Roche_Habenula/raw-data/FASTQ_snRNAseq

cp -r /dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2021-03-15_KMay022421/Br1204_Hb* .
cp -r /dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2021-03-15_KMay022421/Br5558_Hb* .
cp -r /dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2021-10-15_KMay110220/Br1469_Hb* .
cp -r /dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2021-05-21_ASpa050421/3c_k* .
cp -r /dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2021-05-21_ASpa050421/1c_k* .
cp -r /dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2021-05-21_ASpa050421/4c_k* .
cp -r /dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2021-05-21_ASpa050421/2c_k* .
```

These are the available files:

```bash
ls -lh /dcs04/lieber/lcolladotor/pilotHb_LIBD001/Roche_Habenula/raw-data/FASTQ_snRNAseq/*

# /dcs04/lieber/lcolladotor/pilotHb_LIBD001/Roche_Habenula/raw-data/FASTQ_snRNAseq/1c_k_L001_ds.0a7b7047146a4e8580777073634c29ae:
# total 5.9G
# -rwxrwx---+ 1 lcollado lieber_lcolladotor 1.9G Jan 17 20:50 1c-k_S13_L001_R1_001.fastq.gz
# -rwxrwx---+ 1 lcollado lieber_lcolladotor 4.0G Jan 17 20:50 1c-k_S13_L001_R2_001.fastq.gz
# 
# /dcs04/lieber/lcolladotor/pilotHb_LIBD001/Roche_Habenula/raw-data/FASTQ_snRNAseq/1c_k_L002_ds.ed2439f840eb47d08d139e21483725f6:
# total 5.9G
# -rwxrwx---+ 1 lcollado lieber_lcolladotor 1.9G Jan 17 20:52 1c-k_S13_L002_R1_001.fastq.gz
# -rwxrwx---+ 1 lcollado lieber_lcolladotor 4.0G Jan 17 20:52 1c-k_S13_L002_R2_001.fastq.gz
# 
# /dcs04/lieber/lcolladotor/pilotHb_LIBD001/Roche_Habenula/raw-data/FASTQ_snRNAseq/2c_k_L001_ds.01c025c8a9274af583a146419a7b973b:
# total 5.4G
# -rwxrwx---+ 1 lcollado lieber_lcolladotor 1.7G Jan 17 20:57 2c-k_S14_L001_R1_001.fastq.gz
# -rwxrwx---+ 1 lcollado lieber_lcolladotor 3.8G Jan 17 20:58 2c-k_S14_L001_R2_001.fastq.gz
# 
# /dcs04/lieber/lcolladotor/pilotHb_LIBD001/Roche_Habenula/raw-data/FASTQ_snRNAseq/2c_k_L002_ds.c77782c03aef40769707e63d06135d8d:
# total 5.4G
# -rwxrwx---+ 1 lcollado lieber_lcolladotor 1.7G Jan 17 21:00 2c-k_S14_L002_R1_001.fastq.gz
# -rwxrwx---+ 1 lcollado lieber_lcolladotor 3.8G Jan 17 20:59 2c-k_S14_L002_R2_001.fastq.gz
# 
# /dcs04/lieber/lcolladotor/pilotHb_LIBD001/Roche_Habenula/raw-data/FASTQ_snRNAseq/3c_k_L001_ds.19d0c09d1d5b4922add0828bf38422e2:
# total 5.2G
# -rwxrwx---+ 1 lcollado lieber_lcolladotor 1.6G Jan 17 20:46 3c-k_S15_L001_R1_001.fastq.gz
# -rwxrwx---+ 1 lcollado lieber_lcolladotor 3.7G Jan 17 20:46 3c-k_S15_L001_R2_001.fastq.gz
# 
# /dcs04/lieber/lcolladotor/pilotHb_LIBD001/Roche_Habenula/raw-data/FASTQ_snRNAseq/3c_k_L002_ds.ee8551ea791a4794ba2553444b7e88fe:
# total 5.2G
# -rwxrwx---+ 1 lcollado lieber_lcolladotor 1.6G Jan 17 20:48 3c-k_S15_L002_R1_001.fastq.gz
# -rwxrwx---+ 1 lcollado lieber_lcolladotor 3.7G Jan 17 20:48 3c-k_S15_L002_R2_001.fastq.gz
# 
# /dcs04/lieber/lcolladotor/pilotHb_LIBD001/Roche_Habenula/raw-data/FASTQ_snRNAseq/4c_k_L001_ds.1e5dca4b37b6461f9e6a2d460d532323:
# total 5.3G
# -rwxrwx---+ 1 lcollado lieber_lcolladotor 1.7G Jan 17 20:54 4c-k_S16_L001_R1_001.fastq.gz
# -rwxrwx---+ 1 lcollado lieber_lcolladotor 3.6G Jan 17 20:54 4c-k_S16_L001_R2_001.fastq.gz
# 
# /dcs04/lieber/lcolladotor/pilotHb_LIBD001/Roche_Habenula/raw-data/FASTQ_snRNAseq/4c_k_L002_ds.d5822165012f40008c02444152516e33:
# total 5.3G
# -rwxrwx---+ 1 lcollado lieber_lcolladotor 1.7G Jan 17 20:55 4c-k_S16_L002_R1_001.fastq.gz
# -rwxrwx---+ 1 lcollado lieber_lcolladotor 3.6G Jan 17 20:56 4c-k_S16_L002_R2_001.fastq.gz
# 
# /dcs04/lieber/lcolladotor/pilotHb_LIBD001/Roche_Habenula/raw-data/FASTQ_snRNAseq/Br1204_Hb_L001_ds.450f8b4b14ce45d99a81d702fad67d74:
# total 6.4G
# -rwxrwx---+ 1 lcollado lieber_lcolladotor 2.0G Jan 17 18:40 Br1204-Hb_S3_L001_R1_001.fastq.gz
# -rwxrwx---+ 1 lcollado lieber_lcolladotor 4.4G Jan 17 18:42 Br1204-Hb_S3_L001_R2_001.fastq.gz
# 
# /dcs04/lieber/lcolladotor/pilotHb_LIBD001/Roche_Habenula/raw-data/FASTQ_snRNAseq/Br1204_Hb_L002_ds.3ebb12ef36b041faba06599157a53d17:
# total 5.9G
# -rwxrwx---+ 1 lcollado lieber_lcolladotor 1.9G Jan 17 18:42 Br1204-Hb_S3_L002_R1_001.fastq.gz
# -rwxrwx---+ 1 lcollado lieber_lcolladotor 4.1G Jan 17 18:44 Br1204-Hb_S3_L002_R2_001.fastq.gz
# 
# /dcs04/lieber/lcolladotor/pilotHb_LIBD001/Roche_Habenula/raw-data/FASTQ_snRNAseq/Br1469_Hb:
# total 14G
# -rwxrwx---+ 1 lcollado lieber_lcolladotor 2.1G Jan 17 20:42 Br1469_Hb_S13_L001_R1_001.fastq.gz
# -rwxrwx---+ 1 lcollado lieber_lcolladotor 4.6G Jan 17 20:44 Br1469_Hb_S13_L001_R2_001.fastq.gz
# -rwxrwx---+ 1 lcollado lieber_lcolladotor 2.1G Jan 17 20:42 Br1469_Hb_S13_L002_R1_001.fastq.gz
# -rwxrwx---+ 1 lcollado lieber_lcolladotor 4.7G Jan 17 20:43 Br1469_Hb_S13_L002_R2_001.fastq.gz
# 
# /dcs04/lieber/lcolladotor/pilotHb_LIBD001/Roche_Habenula/raw-data/FASTQ_snRNAseq/Br5558_Hb_L001_ds.68bb30be674d4e3da8abd44466ecf985:
# total 7.0G
# -rwxrwx---+ 1 lcollado lieber_lcolladotor 2.2G Jan 17 20:38 Br5558-Hb_S2_L001_R1_001.fastq.gz
# -rwxrwx---+ 1 lcollado lieber_lcolladotor 4.9G Jan 17 20:37 Br5558-Hb_S2_L001_R2_001.fastq.gz
# 
# /dcs04/lieber/lcolladotor/pilotHb_LIBD001/Roche_Habenula/raw-data/FASTQ_snRNAseq/Br5558_Hb_L002_ds.c7f8800daced4fddaee43c22fe006a4c:
# total 6.7G
# -rwxrwx---+ 1 lcollado lieber_lcolladotor 2.1G Jan 17 20:39 Br5558-Hb_S2_L002_R1_001.fastq.gz
# -rwxrwx---+ 1 lcollado lieber_lcolladotor 4.7G Jan 17 20:41 Br5558-Hb_S2_L002_R2_001.fastq.gz
```
