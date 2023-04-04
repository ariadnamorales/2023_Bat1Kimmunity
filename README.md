# 2023_Bat1Kimmunity
Customs scripts used for analyses of 2022-2023 Bat1K immunity project
No new software was developed for this study, thus we provide example commands for analyses or provide links to other sources employed.


## Analyses
### - Genome assembly

#### - Contig assembly
  ```
  ## HiFiasm assemblies
  hifiasm -l2 -o asm_mRhiAff.asm -t 38 *fasta.gz
  awk '/^S/{print ">"$2"\n"$3}' asm_mRhiAff.asm.p_ctg.gfa | fold > asm_mRhiAff.asm.p_ctg.fasta
  awk '/^S/{print ">"$2"\n"$3}' asm_mRhiAff.asm.a_ctg.gfa | fold > asm_mRhiAff.asm.a_ctg.fasta
  
  ## HiCanu assemblies
  canu -p asm_mDorCyc -d hicanu genomeSize=2200m -pacbio-hifi *fasta.gz
  
  ## Canu ONT assemblies
  canu -p asm_mMegSpa_canu2.1_ont -d canu genomeSize=2000m -nanopore *fastq.gz
  
  ## Medake polishing
  medaka_consensus -i merge_ont.fastq.gz -d asm_mMegSpa_canu2.1_ont.contigs.fasta -o polished -t 24 -m r941_prom_hac_g507
  
  ## Purge-dups
  purge_dups/src/split_fa consensus.fasta > split.genome.fasta
  purge_dups/src/calcuts coverage/PB.stat > cutoffs
  minimap2 -I 200G -t 1 -xasm5 -DP split.genome.fasta split.genome.fasta > split.genome.paf
  purge_dups/src/purge_dups -2 -c coverage/PB.base.cov -T cutoffs split.genome.paf > dups.bed
  purge_dups/src/get_seqs -e -p asm_mMegSpa.canu.medaka dups.bed consensus.fasta
  ```

#### - Scaffolding
  ```
  ## Scaff10X
  longranger-2.2.2/longranger mkref contigs.fasta
  longranger-2.2.2/longranger align --reference=refdata-contigs --fastqs=10x_dir --sample=mRhiAff
  scaff10x -nodes 24 -bam possorted_bam.bam contigs.fasta scaff10x.fasta
  longranger-2.2.2/longranger mkref scaff10x.fasta
  longranger-2.2.2/longranger align --reference=refdata-scaff10x --fastqs=10x_dir --sample=mRhiAff
  break10x -nodes 24 -bam possorted_bam.bam scaff10x.fasta break10x.fasta
  
  ## Salsa2
  bwa index genome.fasta
  samtools index genome.fasta
  bwa mem -t24 -B8 genome.fasta hic_R1.fastq.gz | samtools view -@24 -Sb - > hic_R1.bam
  samtools view -h hic_R1.bam | perl filter_five_end.pl | samtools view -@24 - > hic_R1.filtered.bam
  bwa mem -t24 -B8 genome.fasta hic_R2.fastq.gz | samtools view -@24 -Sb - > hic_R2.bam
  samtools view -h hic_R2.bam | perl filter_five_end.pl | samtools view -@24 - > hic_R2.filtered.bam
  perl two_read_bam_combiner.pl hic_R1.filtered.bam hic_R2.filtered.bam | samtools view -@24 -Sb > hic.combined.bam
  bash dedup.sh hic.combined.bam 24
  bedtools bamtobed -i sort_dedup/combined.sort.dp.sort_n.bam > combined.sort.dp.sort_n.bed
  python2.7 salsa2/SALSA-2.2/run_pipeline.py -a genome.fasta -l genome.fasta.fai -e enz.txt -b combined.sort.dp.sort_n.bed -o salsa2_out -m yes -i 50 -p yes
  ```
  
#### - Polishing
  ```
  ## PacBio HiFi polishing
  pbmm2 align mRhiAff.fasta bam.fofn asm_mRhiAff.mapped.bam --sort -j 16 -J 8 --preset CCS -N 1
  singularity exec -B /projects --bind /usr/lib/locale/ deepvariant.sif /opt/deepvariant/bin/run_deepvariant --model_type=PACBIO --ref=mRhiAff.newScaff.fasta --reads=asm_mRhiAff.mapped.bam --output_vcf=mRhiAff.mapped.deepVariant.vcf.gz --num_shards=24
  bcftools view --threads=12 -i 'FILTER=\"PASS\" && GT=\"1/1\"' -o mRhiAff.mapped.deepVariant.HomFiltered.vcf.gz mRhiAff.mapped.deepVariant.vcf.gz
  bcftools index mRhiAff.mapped.deepVariant.HomFiltered.vcf.gz
  bcftools consensus -f mRhiAff.newScaff.fasta -o asm_mRhiAff.deepVariant.fasta mRhiAff.mapped.deepVariant.HomFiltered.vcf.gz
  
  ## 10X polishing
  longranger-2.2.2/longranger mkref genome.fasta
  longranger-2.2.2/longranger align --reference=refdata-genome --fastqs=10x_dir --sample=mRhiAff
  freebayes-1.3.2/freebayes/bin/freebayes -f genome.fasta -g 600 --bam possorted.bam | bcftools view --no-version -Ou > freebayes.bcf
  bcftools view -i -Oz freebayes.bcf > freebayes.vcf
  merfin/build/bin/merfin -polish -sequence genome.fasta -readmers 10x_database.meryl -peak 32 -vcf freebayes.vcf -output merfin.vcf
  bcftools view -Oz merfin.vcf > merfin.vcf.gz && bcftools index merfin.vcf.gz
  bcftools consensus -H 1 -f genome.fasta merfin.vcf.gz > polished_genome.fasta
  ```


#### - Annotation of Transposable Elements
Methods and code as describeb by Osmanski et al.,  [In Press](https://www.biorxiv.org/content/10.1101/2022.12.28.522108v1)

#### - [Repeat masking](http://www.repeatmasker.org) for pairwise genome alignemnts
  ```
  ## build database for RepeatModeler
  BuildDatabase -name $query ${genome_query.fa}
  
  ## run RepeatModeler
  RepeatModeler -threads $cpu -database ${genome_query.fa}
  
  ## run RepeatMasker
  RepeatMasker -pa $cpu -xsmall -lib consensi.fa.classified ${genome_query.fa}
  ```

#### - Pairwise genome alignemnts
  Detailed modified UCSC pipeline [here](https://github.com/hillerlab/GenomeAlignmentTools)

#### - Genome annotation using [TOGA](https://github.com/hillerlab/TOGA)
  ```
  toga.py ${chains.ref.query} ${annotation_ref.bed} ${genome_ref.2bit} ${genome_query.2bit} -i ${isoforms.-reftxt} --cb 3,5 --cjn 500 --u12 ${U12sites_ref.tsv} --ms
  ```
#### - Exon-by-exon codon alignments and cleaning using [extract_codon_alignment.py](https://github.com/hillerlab/TOGA/blob/master/supply/extract_codon_alignment.py) and [HmmCleaner](https://metacpan.org/dist/Bio-MUST-Apps-HmmCleaner/view/bin/HmmCleaner.pl)
   ```
   ## exon-by-exon alignment
   extract_codon_alignments_from_toga.py ${f_listSP} ${annotation_ref.bed} ${trasncriptID} -o ${trasncriptID}_raw.fa -s
   
   ## cleaning
   HmmCleaner.pl -costs ${c1} ${c2} ${c3} ${c4} ${transcriptID}_raw.fa
   ```
  
#### - Phylogenetic and Divergence Time Estimation
  - Gene trees using [raxml](https://cme.h-its.org/exelixis/web/software/raxml/)
  ```
  raxmlHPC-PTHREADS -T ${nThreads} -s ${transcriptID}.fa -m ${sustMod} -N ${reps} -p ${seedSearch} -w ${P_out} -f a -N ${bootstrapReps} --bootstop-perms=${bootstrapReps} -n ${transcriptID} -x ${seedSearch} 
  ```
  
  - Coalescent-based species trees with [ASTRAL](https://github.com/smirarab/ASTRAL)
  ```
  ## concatenate all raxml gene trees
  cat ${P_out_raxml}/RAxML_bestTree.${transcriptID}" >> ${all_raxml_inTrees}
  
  ##run ASTRAL
  ${astral_call} -i ${all_raxml_inTrees} -o ${astral_outTree}
  ```

  - Concatenation-based species trees with [iq-tree](http://www.iqtree.org/)
  ```
  ## concatenate all genes
  ## use https://github.com/ballesterus/Utensils/blob/master/geneStitcher.py
  geneStitcher.py -in ${Path_aln_fasta}/*.fasta
  
  ## run iq-tree
  iqtree -s SuperMatrix.fas -spp Partition.txt -bb 1000  -bnni  -m GTR+I+G -nt 16  -mem 50G
   ```
   
  - Divergence time estimation [treePL](https://github.com/blackrim/treePL)
   ```
   ## use config file1 to optimize parameters
   treePL config.1

   ## see output and generate config.2
   treePL config.2 treefile=tree_rooted.dated.config2.tre
   ```
   

#### - Genome-wide Unbiased Selection Screen using [HYPHY-aBSREL](https://stevenweaver.github.io/hyphy-site/methods/selection-methods/)
  For each transcript, the species tree should be trimmed to kepp only branches represented in alignemnt. We used [tree_doctor](https://github.com/UCSantaCruzComputationalGenomicsLab/phast/blob/master/src/util/tree_doctor.c).
   ```
  ## trimm input tree 
  tree_doctor -a -n -P $(grep '>' ${P_out}/${transcriptID}.fa | sed 's/>//g'| cut -f1 | tr '\n' ',') ${tree} > ${P_out}/${trasncriptID}.prunnedTree.tre

  ## run absrel
  hyphy absrel --alignment ${P_out}/${trasncriptID}.fa --tree ${P_out}/${trasncriptID}.prunnedTree.tre --output ${P_out}/${trasncriptID}.ABSREL.json | tee -a ${P_out}/${trasncriptID}.ABSREL.log
   ```

#### - Gene enrichment analyses using [gprofiler2](https://biit.cs.ut.ee/gprofiler/page/r)
   ```R
   ## Load R and library
   library(gprofiler2)
   
   ## run enrichment per mammalian group
   multi_sig_noBG_tmp_gostres = gost(listGenes_groups, organism = "hsapiens", significant = TRUE)
   
   ## save enrichment summary as tab file
   write.table(sapply(multi_sig_noBG_tmp_gostres$result, FUN = paste), file=out_summaryErich, sep="\t", quote=FALSE, row.names=FALSE)
```

#### -  Correlation between branch length and number of genes under selection

  Negative binomial regression was perfomed using several R libraries, [see code here](https://github.com/lmdavalos/count2branches).

