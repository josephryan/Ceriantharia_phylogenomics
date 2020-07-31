# PLANNED ANALYSES FOR ASSEMBLING THE CERIANTHID PHYLOGENY 
 Principle Investigators: Joseph Ryan and Melissa DeBiasse 
 Draft or Version Number: v.1.0  
 Date: 31 July 2020  
 Note: updates to this document will be tracked through github
 
## 1 INTRODUCTION: BACKGROUND INFORMATION AND SCIENTIFIC RATIONALE  

### 1.1 _Background Information_
Ceriantharians are tube-dwelling anthozoans found in marine reef and benthic communities characterized by ptychocytes, a specialized type of cnidocyte that contains a sticky thread. Cerianthids fire their ptychocytes into the sediment and use the grains of sand and mud that adhere to the thread to build their tube. At the molecular level, cerianthids are characterized by linear, instead of circular, mitochondrial genomes.

### 1.2 _Rationale_ 
The position of Ceriantharia within the broader cnidarian tree of life has long been uncertain and three recent phylogenomic studies (Chang et al. 2015; Zapata et al. 2015; Kayal et al. 2018) had low representation of cerianthid taxa.

### 1.3 _Objectives_   
We will estimate a phylogenetic tree of Ceriantharia using genomic and transcriptomic data from a range of cerianthid and outgroup taxa.

## 2 STUDY DESIGN AND ENDPOINTS  

#### 2.1 The following taxa will be included in the study. Data types are in parentheses.

IN GROUP
Arachnanthus sp. (transcriptome) 
Botruanthus mexicanus (transcriptome) 
Ceriantheopsis americanus (Atlantic) (transcriptome) 
Ceriantheopsis americanus (Gulf of Mexico) (transcriptome) 
Cerianthus borealis (transcriptome) 
Ceriantheomorphe brasiliensis (transcriptome) 
Isarachnanthus maderensis (transcriptome) 
Isarachnanthus nocturnus (transcriptome) 
Pachycerianthus maua (transcriptome) 
Pachycerianthus magnus (draft genome)
OUT GROUP
Nematostella vectensis (transcriptome)
Acropora digitifera (transcriptome)
Protopalythoa variabilis (transcriptome)

#### 2.2 Assemble Pachycerianthus magnus genome

2.2.1 trim Illumina reads

```
java Trimmomatic-0.38/trimmomatic-0.38.jar PE -threads 40 -phred33 R1.fastq.gz R2.fastq.gz P_magn_R1_trim.fq P_magn_R1_unp_trim.fq P_magn_R2_trim.fq P_magn_R2_unp_trim.fq ILLUMINACLIP:/usr/local/Trimmomatic-0.38/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 > trimmo.log 2> trimmo.err &
```

2.2.2 error correct trimmed reads

```
# perform a preliminary assembly (plat -k 45) use it to identify the MT genome
cat P_magn_R1_unp_trim.fq P_magn_R2_unp_trim.fq > P_magn_unp.fq
gzip -9 P_magn_unp.fq
perl /usr/local/allpathslg-44837/src/ErrorCorrectReads.pl PAIRED_READS_A_IN=/P_magn_R1_trim.fq.gz PAIRED_READS_B_IN=P_magn_R2_trim.fq.gz UNPAIRED_READS_IN=P_magn_unp.fq.gz PAIRED_SEP=100 THREADS=46 PHRED_ENCODING=33 READS_OUT=P_magn > ecr.out 2> ecr.err &
```

2.2.3 remove mitochondrial reads

```
FastqSifter --out=_P_mag_sifted --fasta=P_magnus_mtdna.fa --right=P_magn.paired.B.fastq --left=P_magn.paired.A.fastq --unp=P_magn.unpaired.fastq --threads=40 --savereads > fqs.out 2> fqs.err &
```

2.2.4 assemble genomes using multiple kmer values

```
plat.pl --out=P_magn_31 --k=31 --m=450 --left=P_mag_sifted.filtered.A.fq --right=P_mag_sifted.filtered.B.fq --unp=P_mag_sifted.filtered.unp.fq --threads=50 > plat.out 2> plat.err &
plat.pl --out=P_magn_53 --k=53 --m=450 --left=P_mag_sifted.filtered.A.fq --right=P_mag_sifted.filtered.B.fq --unp=P_mag_sifted.filtered.unp.fq --threads=50 > plat.out 2> plat.err &
plat.pl --out=P_magn_69 --k=69 --m=450 --left=P_mag_sifted.filtered.A.fq --right=P_mag_sifted.filtered.B.fq --unp=P_mag_sifted.filtered.unp.fq --threads=50 > plat.out 2> plat.err &
plat.pl --out=P_magn_73 --k=73 --m=450 --left=P_mag_sifted.filtered.A.fq --right=P_mag_sifted.filtered.B.fq --unp=P_mag_sifted.filtered.unp.fq --threads=50 > plat.out 2> plat.err &
```

2.2.5 concatenate genomes

```
cat P_magn_31_out_gapClosed.fa P_magn_53_out_gapClosed.fa P_magn_69_out_gapClosed.fa P_magn_73_out_gapClosed.fa > P_magn_genomes.fa 
```

#### 2.2 Translate the concatenated nucleotide genome sequence into amino acids in TransDecoder v3.0.1. We set the –m flag to 50 and used the results from BLAST searches to inform the final TransDecoder prediction step

```
TransDecoder.LongOrfs -t P_magn_genomes.fa -m 50 > td.out 2> td.err
```

```
blastp -query longest_orfs.pep -db swissprot -max_target_seqs 1 -outfmt 6 -evalue 1e-5 -num_threads [#] > blastp.out 2> blastp.err 
```

```
TransDecoder.Predict -t P_magn_genomes.fa --retain_blastp_hits out.blastp.out > tdp.out 2> tdp.err
```

#### 2.3 Build single copy gene matrix for cerianthid and outgroup taxa 

2.3.1 identify orthologs among taxa using OrthoFinder v2.2.3 

```
orthofinder -f [dir_w_protein_fasta_files] -S diamond -M msa -os > of.out 2> of.err &

```

2.3.2 select orthogroups that have P. magnus and at least 6 other taxa
```
pmag_orthogroup_filter.pl
```

2.3.3 align orthogroup sequences and estimate a gene tree for each orthogroup
```
mafft [OG.fa] > OG.mafft
iqtree-omp -s OG.mafft -pre OG -nt AUTO -m TEST -bb 1000 > OG.iq.stdout 2> OG.iq.err
```

2.3.4 Prune each gene tree into maximally inclusive subtrees with no more than one sequence per taxon. Each of these subtrees is considered as a set of putative orthologs. (method used by Agalma, Dunn et al. 2013) We will do manually as it will be a small number of trees.


2.3.5 concatenate the single-copy loci filtered from step 2.3.4 to create a matrix and partition file for use in downstream phylogenomic analyses using ```fasta2phylomatrix``` (available in https://github.com/josephryan/RyanLabPhylogenomicTools).

```
fasta2phylomatrix  {--fasta=FILE1 --fasta=FILE2 ... or --dir=DIR_OF_FASTAFILES} [--nexpartition=OUT_PARTITION_FILE] [--raxpartition=OUT_PARTITION_FILE] > cerianthid_matrix
```

#### 2.4 Estimate cerianthid species phylogenies using maximum likelihood

```
iqtree-omp -s cerianthid_matrix.fa -pre cerianthid -spp cerianthid.nex -nt AUTO -m TEST -bb 1000 > iq.stdout 2> iq.err
``` 

#### 2.5 Compare topologies of alternative matrices to evaluate variation in phylogenetic signal among gene sets

2.5.1 estimate phylogeny for cerianthid taxa using 748 gene matrix from DeBiasse et al., a study to infer the cnidarian tree of life (for details, see https://github.com/josephryan/DeBiasse_cnidophylogenomics). This matrix does not include P. magnus

```
iqtree-omp -s cerianthid_748_matrix.fa -pre cerianthid_748 -spp cerianthid_748.nex -nt AUTO -m TEST -bb 1000 > iq.stdout 2> iq.err
```

2.5.2 compare the cerianthid relationships in the cnidarian species tree generated in DeBiasse et al. to the relationships in the cerianthid species tree generated in step 2.5.1. These trees have the same set of cerianthid taxa and gene set, but  differ in outgroup taxa and number of outgroups. We expect the relationships among cerianthids in these trees to be the same, indicating that the choice of outgroups has no effect on the topology of the tree. If the relationships among cerianthids in the trees differs, we will update the phylotocol with a strategy to further investigate the discordance.

2.5.3 compare the cerianthid relationships in the cerianthid_748 species tree with the cerianthid relationships in the tree generated in 2.4 (after pruning P. magnus). These trees have the same taxa, but a different gene set. If the topologies of these trees differ, we will conclude that the smaller gene matrix we generated in this study lacks phylogenetic signal and therefore, we will not be able to use it to place P. magnus in the cerianthid tree.


#### 3 Work completed to-date
July 21 2020 steps 2.2 and 2.3.1 - 2.3.3 have been completed


## 4 PROGRAMS REFERENCED  

Dunn Dunn CW, Howison M, Zapata F. Agalma: an automated phylogenomics workflow. BMC bioinformatics. 2013 Dec;14(1):1-9.

Emms, D. M., & Kelly, S. (2015). OrthoFinder: solving fundamental biases in whole genome comparisons dramatically improves orthogroup inference accuracy. Genome Biology, 16(1), 157. 

Nguyen, L. T., Schmidt, H. A., von Haeseler, A., & Minh, B. Q. (2014). IQ-TREE: a fast and effective stochastic algorithm for estimating maximum-likelihood phylogenies. Molecular Biology and Evolution, 32(1), 268-274.

TransDecoder: https://transdecoder.github.io/

Yamada, K. D., Tomii, K., & Katoh, K. (2016). Application of the MAFFT sequence alignment program to large data—reexamination of the usefulness of chained guide trees. Bioinformatics, 32(21), 3246-3251.

## APPENDIX

Version : Date : Significant Revisions  
1.1
1.2  
1.3  
1.4 
